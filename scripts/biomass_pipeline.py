#!/usr/bin/env python
"""End-to-end BIOMASS demo: download -> GSLC -> interferograms.

Ties the pieces together for a repeat-pass stack over one track/frame:

1. download N L1a SCS scenes from ESA-MAAP (``download_maap`` helpers),
2. fetch a Copernicus DEM for the footprint (``sardem``) and reproject to UTM,
3. geocode every scene's chosen polarization to a *single shared geogrid*
   (so the GSLCs are coregistered by construction) and save them as GeoTIFFs,
4. form sequential interferograms + coherence and render a PNG.

The saved GSLCs (``<work>/gslc/biomass_<pol>_<date>.tif``) are directly
consumable by dolphin for phase linking / unwrapping.

Usage
-----
    python scripts/biomass_pipeline.py \
        --bbox -66.5 2.5 -64.5 4.5 --start 2026-04-22 --end 2026-04-30 \
        --track T007 --frame F004 --max 3 --pol HH \
        --work /path/to/work --spacing 30 --looks 5
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import download_maap  # noqa: E402


def _utm_epsg(lon: float, lat: float) -> int:
    zone = int((lon + 180) // 6) + 1
    return (32600 if lat >= 0 else 32700) + zone


def fetch_dem(bbox: list[float], work: pathlib.Path) -> pathlib.Path:
    """Download a Copernicus DEM for the bbox and reproject to UTM (30 m)."""
    w, s, e, n = bbox
    dem_ll = work / "dem_ll.tif"
    dem_utm = work / "dem_utm.tif"
    if not dem_utm.exists():
        subprocess.run(
            ["sardem", "--bbox", str(w), str(s), str(e), str(n),
             "--data-source", "COP", "-o", str(dem_ll)],
            check=True,
        )
        epsg = _utm_epsg((w + e) / 2, (s + n) / 2)
        subprocess.run(
            ["gdalwarp", "-t_srs", f"EPSG:{epsg}", "-tr", "30", "30",
             "-r", "bilinear", "-overwrite", str(dem_ll), str(dem_utm)],
            check=True,
        )
    return dem_utm


def geocode_stack(
    product_dirs: list[pathlib.Path],
    dem_utm: pathlib.Path,
    pol: str,
    spacing: float,
    gslc_dir: pathlib.Path,
) -> list[pathlib.Path]:
    """Geocode every scene to one shared geogrid; save GSLC GeoTIFFs."""
    import isce3
    from osgeo import gdal, osr

    from biomass_reader import BiomassSlc

    gslc_dir.mkdir(parents=True, exist_ok=True)
    dem = isce3.io.Raster(str(dem_utm))
    epsg = dem.get_epsg()
    ellipsoid = isce3.core.make_projection(epsg).ellipsoid

    scenes = sorted(
        (BiomassSlc.from_dir(d, pol) for d in product_dirs),
        key=lambda s: s.sensing_start,
    )
    ref = scenes[0]
    geogrid = isce3.product.bbox_to_geogrid(
        ref.radar_grid, ref.orbit, ref.doppler, spacing, -spacing, epsg
    )
    print(f"shared geogrid {geogrid.length} x {geogrid.width}  EPSG:{epsg}")

    out_paths = []
    for slc in scenes:
        geo = np.full(
            (geogrid.length, geogrid.width), np.nan + 1j * np.nan, dtype=np.complex64
        )
        t = time.time()
        isce3.geocode.geocode_slc(
            geo, slc.read_complex(), dem, slc.radar_grid, geogrid, slc.orbit,
            slc.doppler, isce3.core.LUT2d(), ellipsoid, 1.0e-8, 25, flatten=True,
        )
        date = slc.sensing_start.strftime("%Y%m%d")
        out = gslc_dir / f"biomass_{pol}_{date}.tif"
        ds = gdal.GetDriverByName("GTiff").Create(
            str(out), geogrid.width, geogrid.length, 1, gdal.GDT_CFloat32,
            options=["COMPRESS=DEFLATE", "TILED=YES"],
        )
        gg = geogrid
        ds.SetGeoTransform(
            (gg.start_x, gg.spacing_x, 0, gg.start_y, 0, gg.spacing_y)
        )
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(np.nan_to_num(geo))
        ds.FlushCache()
        ds = None
        print(f"  geocoded {date} in {time.time() - t:.1f}s -> {out.name}")
        out_paths.append(out)
    return out_paths


def make_interferograms(
    gslc_paths: list[pathlib.Path], looks: int, out_png: pathlib.Path
):
    """Form sequential interferograms + coherence and render a PNG."""
    import matplotlib
    import rasterio

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402
    from sarlet.plot import render_phase
    from sarlet.utils import take_looks

    def read(p):
        with rasterio.open(p) as src:
            return src.read(1)

    arrs = [read(p) for p in gslc_paths]
    dates = [p.stem.split("_")[-1] for p in gslc_paths]
    n_pairs = len(arrs) - 1
    _, ax = plt.subplots(n_pairs, 3, figsize=(15, 6 * n_pairs), squeeze=False)
    for i in range(n_pairs):
        g1, g2 = arrs[i], arrs[i + 1]
        looked = take_looks(g1 * g2.conj(), looks, looks, "nansum")
        p1 = take_looks(np.abs(g1) ** 2, looks, looks, "nansum")
        p2 = take_looks(np.abs(g2) ** 2, looks, looks, "nansum")
        coh = np.clip(np.nan_to_num(np.abs(looked) / np.sqrt(p1 * p2 + 1e-12)), 0, 1)
        mask = (p1 > 0) & (p2 > 0)
        amp = 20 * np.log10(np.abs(looked) + 1e-6)
        amp[~mask] = np.nan
        vmn, vmx = np.nanpercentile(amp, [5, 95])
        ax[i, 0].imshow(amp, cmap="gray", vmin=vmn, vmax=vmx)
        ax[i, 0].set_title(f"{dates[i]} x {dates[i+1]}  amplitude")
        ax[i, 1].imshow(render_phase(looked, coherence=coh))
        ax[i, 1].set_title("wrapped phase (coh-modulated)")
        c = coh.copy()
        c[~mask] = np.nan
        im = ax[i, 2].imshow(c, cmap="magma", vmin=0, vmax=1)
        ax[i, 2].set_title(f"coherence (mean {np.nanmean(c):.2f})")
        for a in ax[i]:
            a.axis("off")
        plt.colorbar(im, ax=ax[i, 2], shrink=0.7)
    plt.suptitle("BIOMASS P-band repeat-pass interferograms (geocoded)", fontsize=14)
    plt.tight_layout()
    plt.savefig(out_png, dpi=95, bbox_inches="tight")
    print(f"wrote {out_png}")


def main() -> None:  # noqa: D103
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", nargs=4, type=float, required=True,
                   metavar=("W", "S", "E", "N"))
    p.add_argument("--start", required=True)
    p.add_argument("--end", required=True)
    p.add_argument("--track", default=None)
    p.add_argument("--frame", default=None)
    p.add_argument("--max", type=int, default=3)
    p.add_argument("--pol", default="HH")
    p.add_argument("--spacing", type=float, default=30.0)
    p.add_argument("--looks", type=int, default=5)
    p.add_argument("--work", type=pathlib.Path, required=True)
    p.add_argument("--credentials", type=pathlib.Path,
                   default=pathlib.Path("~/.maap/credentials.txt"))
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)

    # 1. download
    creds = download_maap.read_credentials(args.credentials)
    token = download_maap.get_access_token(creds)
    feats = download_maap.search(args.bbox, args.start, args.end)
    if args.track:
        feats = [f for f in feats if f"_{args.track}_" in f["id"]]
    if args.frame:
        feats = [f for f in feats if f"_{args.frame}_" in f["id"]]
    feats.sort(key=lambda f: f["properties"]["datetime"])
    feats = feats[: args.max]
    product_dirs = []
    for f in feats:
        zip_path = args.work / f"{f['id']}.zip"
        out_dir = args.work / f["id"]
        if not out_dir.exists():
            download_maap.download(f["assets"]["product"]["href"], token, zip_path)
            import zipfile
            with zipfile.ZipFile(zip_path) as z:
                z.extractall(args.work)
        product_dirs.append(out_dir)

    # 2-4
    dem_utm = fetch_dem(args.bbox, args.work)
    gslc_paths = geocode_stack(product_dirs, dem_utm, args.pol, args.spacing,
                               args.work / "gslc")
    make_interferograms(gslc_paths, args.looks, args.work / "interferograms.png")
    print("\nGSLCs ready for dolphin:")
    for g in gslc_paths:
        print("  ", g)


if __name__ == "__main__":
    main()
