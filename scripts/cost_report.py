#!/usr/bin/env python3
"""Runtime and CPU cost per pipeline step from a Nextflow trace, extrapolated to larger runs.

Usage (after a run; needs pandas, e.g. inside the analysis image):
  python scripts/cost_report.py --trace results/pilot/pipeline_info/trace_*.tsv --genomes 2000 \
      --scale 20000 54963 --out reports/pilot/cost.md

Per-genome steps (fetch, barrnap, genes/16S, eggNOG) are extrapolated linearly in genome count.
Whole-set steps are listed as measured: their scaling differs (alignment ~linear, Aim 1 models
~n log n to n^2, Aim 2 ~n^2 per test) and is stated next to each.
"""
import argparse
import re

import pandas as pd

PER_GENOME = {"FETCH_GENOMES", "BARRNAP", "PROCESS_GENOMES", "EGGNOG"}
SCALING = {"ALIGN_16S": "~linear", "BUILD_FEATURES": "~linear", "BUILD_TARGETS": "~linear",
           "AIM1_PREDICT": "~n log n (RF) to n^2 (k-NN)", "AIM2_PYSEER": "~n^2 per test (LMM)",
           "PICRUST2": "~linear + fixed reference cost", "PATRISTIC": "n^2 memory (D063)"}


def seconds(s):
    """Nextflow duration string ('1h 2m 3s', '850ms', '4.5s') -> seconds."""
    if not isinstance(s, str) or s == "-":
        return 0.0
    units = {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 0.001}
    return sum(float(v) * units[u] for v, u in re.findall(r"([\d.]+)(ms|d|h|m|s)", s))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace", nargs="+", required=True)
    ap.add_argument("--genomes", type=int, required=True)
    ap.add_argument("--scale", type=int, nargs="+", default=[20000, 54963])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    t = pd.concat(pd.read_csv(f, sep="\t") for f in a.trace)
    t = t[t.status.isin(["COMPLETED", "CACHED"])].drop_duplicates("hash", keep="last")
    t["step"] = t["name"].str.split(" ").str[0].str.split(":").str[-1]
    t["wall_s"] = t["realtime"].map(seconds)
    t["cpu_h"] = t["wall_s"] * pd.to_numeric(t["%cpu"].astype(str).str.rstrip("%"), errors="coerce").fillna(100) / 100 / 3600
    s = t.groupby("step").agg(tasks=("hash", "size"), wall_h=("wall_s", lambda x: x.sum() / 3600),
                              cpu_h=("cpu_h", "sum")).reset_index()
    s["scaling"] = s.step.map(lambda x: "linear (per genome)" if x in PER_GENOME else SCALING.get(x, "small"))
    for n in a.scale:
        s[f"cpu_h_at_{n}"] = s.apply(lambda r: r.cpu_h * n / a.genomes if r.step in PER_GENOME else float("nan"), axis=1)
    lines = ["# Cost report", "", f"Measured on {a.genomes} genomes. Wall hours are summed task time, "
             "not elapsed time (tasks overlap).", "",
             "| " + " | ".join(s.columns) + " |", "|" + "---|" * s.shape[1]]
    lines += ["| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |"
              for r in s.itertuples(index=False)]
    pg = s[s.step.isin(PER_GENOME)]
    lines += ["", f"Per-genome steps: {pg.cpu_h.sum() / a.genomes * 3600:.1f} CPU-seconds per genome."]
    for n in a.scale:
        cpu = pg.cpu_h.sum() * n / a.genomes
        lines.append(f"- {n} genomes: ~{cpu:.0f} CPU-hours for per-genome steps (~{cpu / 32:.0f} h on 32 cores).")
    open(a.out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
