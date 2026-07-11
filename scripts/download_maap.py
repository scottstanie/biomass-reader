#!/usr/bin/env python
"""Download BIOMASS L1a SCS products from the ESA-MAAP STAC catalog.

Catalog *search* is public; *download* needs a MAAP Bearer token, obtained by
exchanging a personal offline token. Provide credentials in a KEY=VALUE file
(default ``~/.maap/credentials.txt``) with::

    CLIENT_ID=offline-token
    CLIENT_SECRET=...
    OFFLINE_TOKEN=<offline token from portal.maap.eo.esa.int>

Usage
-----
    python scripts/download_maap.py \
        --bbox -66.5 2.5 -64.5 4.5 \
        --start 2026-04-22 --end 2026-04-30 \
        --frame F004 --track T007 \
        --dest /path/to/dir --max 3 --unzip
"""

from __future__ import annotations

import argparse
import pathlib
import time
import zipfile

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

IAM = "https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token"
CATALOG = "https://catalog.maap.eo.esa.int/catalogue"


def _session() -> requests.Session:
    """Return a requests session that retries connection resets with backoff.

    The MAAP endpoints intermittently reset the TLS connection; urllib3's
    ``Retry`` transparently re-establishes it for connect/read errors and 5xx
    responses on both POST (token/search) and GET (download).
    """
    retry = Retry(
        total=6,
        connect=6,
        read=6,
        backoff_factor=0.5,
        backoff_max=20,
        status_forcelist=(500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        raise_on_status=False,
    )
    s = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    return s


SESSION = _session()


def read_credentials(path: pathlib.Path) -> dict[str, str]:
    """Parse a KEY=VALUE MAAP credentials file."""
    kv = {}
    for line in path.expanduser().read_text().splitlines():
        line = line.strip()
        if line and "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()
    for req in ("CLIENT_ID", "CLIENT_SECRET", "OFFLINE_TOKEN"):
        assert req in kv, f"{req} missing from credentials file {path}"
    return kv


def get_access_token(creds: dict[str, str], retries: int = 4) -> str:
    """Exchange the offline token for a short-lived access token."""
    data = {
        "client_id": creds["CLIENT_ID"],
        "client_secret": creds["CLIENT_SECRET"],
        "grant_type": "refresh_token",
        "refresh_token": creds["OFFLINE_TOKEN"],
    }
    for attempt in range(retries):
        try:
            r = SESSION.post(IAM, data=data, timeout=30)
            r.raise_for_status()
            return r.json()["access_token"]
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            print(f"  token exchange retry {attempt + 1}: {e}")
            time.sleep(2**attempt)
    raise RuntimeError("unreachable")


def search(
    bbox: list[float], start: str, end: str, limit: int = 200
) -> list[dict]:
    """Search the BiomassLevel1a SCS collection (public, no auth)."""
    body = {
        "collections": ["BiomassLevel1a"],
        "filter-lang": "cql2-json",
        "filter": {"op": "=", "args": [{"property": "productType"}, "S1_SCS__1S"]},
        "bbox": bbox,
        "datetime": f"{start}T00:00:00Z/{end}T00:00:00Z",
        "limit": limit,
    }
    r = SESSION.post(f"{CATALOG}/search", json=body, timeout=60)
    r.raise_for_status()
    return r.json()["features"]


def download(
    url: str, token: str, dest: pathlib.Path, retries: int = 4
) -> pathlib.Path:
    """Stream a product zip to ``dest`` with simple resume-on-retry."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        headers = {"Authorization": f"Bearer {token}"}
        pos = dest.stat().st_size if dest.exists() else 0
        if pos:
            headers["Range"] = f"bytes={pos}-"
        mode = "ab" if pos else "wb"
        try:
            with SESSION.get(url, headers=headers, stream=True, timeout=120) as r:
                if r.status_code not in (200, 206):
                    r.raise_for_status()
                total = int(r.headers.get("content-length", 0)) + pos
                with open(dest, mode) as f:
                    done = pos
                    for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                        done += f.write(chunk)
                        if total:
                            pct = 100 * done / total
                            print(
                                f"\r  {dest.name[:48]}  {done/1e6:7.1f}"
                                f"/{total/1e6:.1f} MB ({pct:4.1f}%)",
                                end="",
                            )
            print()
            return dest
        except requests.RequestException as e:
            print(f"\n  download retry {attempt + 1}: {e}")
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {url}")


def main() -> None:  # noqa: D103
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", nargs=4, type=float, required=True,
                   metavar=("W", "S", "E", "N"))
    p.add_argument("--start", required=True)
    p.add_argument("--end", required=True)
    p.add_argument("--track", default=None, help="e.g. T007")
    p.add_argument("--frame", default=None, help="e.g. F004")
    p.add_argument("--dest", type=pathlib.Path, required=True)
    p.add_argument("--max", type=int, default=3)
    p.add_argument("--credentials", type=pathlib.Path,
                   default=pathlib.Path("~/.maap/credentials.txt"))
    p.add_argument("--access-token-file", type=pathlib.Path, default=None,
                   help="use a pre-fetched access token instead of exchanging "
                        "the offline token (avoids the IAM host)")
    p.add_argument("--unzip", action="store_true")
    args = p.parse_args()

    if args.access_token_file:
        token = args.access_token_file.expanduser().read_text().strip()
        print("using pre-fetched access token")
    else:
        creds = read_credentials(args.credentials)
        token = get_access_token(creds)
        print("access token acquired")

    for attempt in range(6):
        try:
            feats = search(args.bbox, args.start, args.end)
            break
        except requests.RequestException as e:
            print(f"  search retry {attempt + 1}: {e}")
            time.sleep(min(2**attempt, 20))
    else:
        raise RuntimeError("search failed after retries")
    if args.track:
        feats = [f for f in feats if f"_{args.track}_" in f["id"]]
    if args.frame:
        feats = [f for f in feats if f"_{args.frame}_" in f["id"]]
    feats.sort(key=lambda f: f["properties"]["datetime"])
    feats = feats[: args.max]
    print(f"{len(feats)} products selected:")
    for f in feats:
        print("  ", f["id"])

    args.dest.mkdir(parents=True, exist_ok=True)
    for f in feats:
        url = f["assets"]["product"]["href"]
        zip_path = args.dest / f"{f['id']}.zip"
        if zip_path.exists() and zipfile.is_zipfile(zip_path):
            print(f"  exists, skipping: {zip_path.name}")
        else:
            download(url, token, zip_path)
        if args.unzip:
            out_dir = args.dest / f["id"]
            if not out_dir.exists():
                print(f"  unzipping {zip_path.name}")
                with zipfile.ZipFile(zip_path) as z:
                    z.extractall(args.dest)
    print("done")


if __name__ == "__main__":
    main()
