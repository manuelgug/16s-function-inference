#!/usr/bin/env python3
"""Download genome assemblies from the NCBI FTP tree over HTTPS, verifying NCBI's md5 checksums.

Usage: fetch_genomes.py GENOMES_TSV OUTDIR [--threads 4]
Writes OUTDIR/<gtdb accession>.fna.gz per genome and fetch_status.tsv (accession, status, detail).
A failed genome is reported, never fatal: downstream steps use what was fetched.
If an RS_ (GCF) assembly is gone, the paired GenBank (GCA) assembly is tried.
"""
import argparse
import csv
import hashlib
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
UA = {"User-Agent": "16s-function-inference/1.0 (+https://github.com/manuelgug/16s-function-inference)"}


def get(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        time.sleep(2 ** i * 3)
    raise RuntimeError(f"{url}: {err}")


def assembly_dir(acc):
    """GCF_000005845.2 -> URL of its directory, or None if NCBI has no such assembly."""
    m = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+", acc)
    if not m:
        raise ValueError(f"not an NCBI assembly accession: {acc}")
    parent = f"{BASE}/{m[1]}/{m[2]}/{m[3]}/{m[4]}/"
    try:
        listing = get(parent).decode()
    except urllib.error.HTTPError:
        return None
    names = sorted(set(re.findall(rf'href="({re.escape(acc)}_[^"/]+)/"', listing)))
    return parent + names[0] + "/" if names else None


def fetch_one(gtdb_acc, ncbi_acc, outdir):
    dest = outdir / f"{gtdb_acc}.fna.gz"
    candidates = [ncbi_acc] + ([ncbi_acc.replace("GCF_", "GCA_")] if ncbi_acc.startswith("GCF_") else [])
    for acc in candidates:
        d = assembly_dir(acc)
        if d is None:
            continue
        name = d.rstrip("/").rsplit("/", 1)[1] + "_genomic.fna.gz"
        md5s = dict(
            (line.split()[1].lstrip("./"), line.split()[0])
            for line in get(d + "md5checksums.txt").decode().splitlines() if line.strip()
        )
        data = get(d + name)
        if name in md5s and hashlib.md5(data).hexdigest() != md5s[name]:
            return gtdb_acc, "failed", f"md5 mismatch for {name}"
        tmp = dest.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.rename(dest)
        return gtdb_acc, "ok", acc
    return gtdb_acc, "failed", f"no NCBI assembly directory for {' or '.join(candidates)}"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("genomes_tsv")
    ap.add_argument("outdir")
    ap.add_argument("--threads", type=int, default=4, help="parallel downloads; keep low to be polite to NCBI")
    a = ap.parse_args(argv)
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(a.genomes_tsv) as f:
        rows = [(r["accession"], r["ncbi_accession"]) for r in csv.DictReader(f, delimiter="\t")]

    def safe(row):
        try:
            return fetch_one(*row, outdir)
        except Exception as e:  # report and continue; one bad genome must not sink the batch
            return row[0], "failed", str(e)[:300]

    with ThreadPoolExecutor(a.threads) as ex:
        results = list(ex.map(safe, rows))
    with open("fetch_status.tsv", "w") as f:
        f.write("accession\tstatus\tdetail\n")
        for r in results:
            f.write("\t".join(r) + "\n")
    n_ok = sum(r[1] == "ok" for r in results)
    print(f"fetched {n_ok}/{len(rows)} genomes", file=sys.stderr)


if __name__ == "__main__":
    main()
