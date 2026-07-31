"""Download NHANES XPT components with retries and atomic replacement."""
from __future__ import annotations

import csv
import hashlib
import os
import time
from pathlib import Path

import requests

ROOT = Path.cwd() / "work"
RAW = ROOT / "data" / "raw"

FILES = {
    "pre": {
        "DEMO": "P_DEMO", "DR1": "P_DR1TOT", "DR2": "P_DR2TOT", "LUX": "P_LUX",
        "BMX": "P_BMX", "BPX": "P_BPXO", "BPQ": "P_BPQ", "GHB": "P_GHB",
        "BIO": "P_BIOPRO", "HDL": "P_HDL", "DIQ": "P_DIQ", "ALQ": "P_ALQ",
        "SMQ": "P_SMQ", "PAQ": "P_PAQ", "HEPBD": "P_HEPBD",
        "HEPC": "P_HEPC", "HEQ": "P_HEQ",
    },
    "post": {
        "DEMO": "DEMO_L", "DR1": "DR1TOT_L", "DR2": "DR2TOT_L", "LUX": "LUX_L",
        "BMX": "BMX_L", "BPX": "BPXO_L", "BPQ": "BPQ_L", "GHB": "GHB_L",
        "BIO": "BIOPRO_L", "HDL": "HDL_L", "DIQ": "DIQ_L", "ALQ": "ALQ_L",
        "SMQ": "SMQ_L", "PAQ": "PAQ_L", "HEPBD": "HEPBD_L",
        "HEPC": "HEPC_L", "HEQ": "HEQ_L",
    },
}

BASES = {
    "pre": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles",
    "post": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles",
}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def valid_xpt(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 10_000:
        return False
    with path.open("rb") as f:
        head = f.read(80)
    return b"HEADER RECORD" in head


def download(url: str, dest: Path) -> None:
    if valid_xpt(dest):
        return
    tmp = dest.with_suffix(".part")
    for attempt in range(1, 6):
        try:
            with requests.get(url, stream=True, timeout=(20, 120)) as response:
                response.raise_for_status()
                with tmp.open("wb") as f:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            f.write(chunk)
            if not valid_xpt(tmp):
                raise RuntimeError(f"Downloaded file failed XPT header/size check: {tmp}")
            os.replace(tmp, dest)
            return
        except Exception:
            if tmp.exists():
                tmp.unlink()
            if attempt == 5:
                raise
            time.sleep(attempt * 2)


rows = []
for period, components in FILES.items():
    folder = RAW / period
    folder.mkdir(parents=True, exist_ok=True)
    for key, stem in components.items():
        destination = folder / f"{stem}.xpt"
        url = f"{BASES[period]}/{stem}.xpt"
        print(f"{period}:{key} {stem}")
        download(url, destination)
        rows.append({
            "period": period,
            "key": key,
            "stem": stem,
            "url": url,
            "bytes": destination.stat().st_size,
            "md5": md5(destination),
        })

with (RAW / "download_manifest.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Downloaded/verified {len(rows)} files.")
