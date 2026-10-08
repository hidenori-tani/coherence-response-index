"""Shared data-acquisition helpers: disk guard + md5-verified ENA download.

Streaming policy (disk ~24-29 GB free; mRNA layer is ~48 GB total): never keep
more than one run's fastq on disk. Each run is downloaded, quantified, then its
fastq is deleted before the next run.
"""
import hashlib
import shutil
import subprocess
from pathlib import Path

from tmc.config import load_accessions

ACC = "config/accessions.tsv"
MIN_FREE_GB = 6.0  # abort before a download if free space drops below this


def free_gb(path="."):
    return shutil.disk_usage(path).free / 1e9


def assert_space(path="."):
    f = free_gb(path)
    if f < MIN_FREE_GB:
        raise SystemExit(f"ABORT: only {f:.1f} GB free (< {MIN_FREE_GB} GB guard). "
                         "Free space and rerun; partial results are kept.")
    return f


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest, want_md5=None, attempts=4):
    """Download via curl (system cert store; macOS Python urllib rejects the ENA
    cert chain). Retries on transient truncation: re-download if md5 mismatches."""
    dest = Path(dest)
    if dest.exists() and want_md5 and md5(dest) == want_md5:
        return dest  # cached + verified
    full = url if url.startswith("http") else "https://" + url
    last = None
    for i in range(attempts):
        try:
            subprocess.run(["curl", "-fsSL", "--retry", "5", "--retry-all-errors",
                            "--retry-delay", "3", "-o", str(dest), full], check=True)
        except subprocess.CalledProcessError as e:
            last = f"curl exit {e.returncode}"
            dest.unlink(missing_ok=True)
            print(f"  download retry {i+1}/{attempts} for {dest.name} ({last})")
            continue
        if not want_md5:
            return dest
        got = md5(dest)
        if got == want_md5:
            return dest
        last = f"md5 {got[:8]}"
        dest.unlink(missing_ok=True)
        print(f"  md5 retry {i+1}/{attempts} for {dest.name} (got {got[:8]}...)")
    raise SystemExit(f"download failed {dest.name} after {attempts} tries: {last} != {want_md5}")


def layer_rows(layer):
    return [r for r in load_accessions(ACC) if r["layer"] == layer]
