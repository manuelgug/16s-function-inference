#!/usr/bin/env python3
"""Benchmark evaluation: every method on the same held-out genomes, the pre-registered win rule,
and the "when and why" analyses. Writes tables, figures and the Markdown + HTML report.

Usage: bench_eval.py --data bench_data/ --splits splits.tsv --tasks task_dir [...] --phase dev|confirm
                     [--thresholds picrust2_thresholds.tsv] [--stacker-importance stacker_importance.tsv]
                     [--trace trace.tsv ...] [--boot 2000] [--seed 42] [--cpus 16] --outdir report/

Each task dir (one level x fold x arm) holds task.tsv, novelty.tsv, nsti.tsv and pred_<method>.npz.
Co-primary endpoints (BENCHMARK.md), full-length 16S, levels genus..class, s16fun minus PICRUSt2's best setting,
95 % CIs by resampling genera: per-KO AUROC superiority (>= +0.05, CI > 0) and KO F1 non-inferiority
(CI lower bound >= -0.01). Win = both, at every level.
"""
import argparse
from pathlib import Path

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from bench_figures import (METHODS, PRIMARY_LEVELS, fig_arms, fig_methods_diagram, fig_novelty,  # noqa: E402
                           fig_ordination, fig_per_ko, fig_phyla, fig_primary, fig_resources)
from bench_ordination import ordination  # noqa: E402
from benchlib import WeightedAUROC, cluster_bootstrap, ko_scores, load_data, load_pred  # noqa: E402
from build_targets import ModuleScorer  # noqa: E402

LEVELS = ["random", "genus", "family", "order", "class"]
P2_METHODS = ["mp", "emp_prob", "pic", "subtree_average"]
THRESHOLD_GRID = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5]
AUROC_MARGIN = 0.05        # co-primary (a): per-KO AUROC superiority margin
F1_NONINFERIORITY = 0.01   # co-primary (b): KO F1 non-inferiority margin


# ---------- per task ----------

def task_metrics(task, truth_Y, truth_mod, kos, thresholds, scorer):
    """Per-genome metrics for every method of one task."""
    info = pd.read_csv(task / "task.tsv", sep="\t")
    if info.empty:
        return pd.DataFrame()
    level, fold, arm, phase = info.iloc[0][["level", "fold", "arm", "phase"]]
    q = info.accession.tolist()
    Y = truth_Y.loc[q].to_numpy()
    M_true = truth_mod.loc[q]
    rows = []
    preds = {f.stem.removeprefix("pred_"): f for f in sorted(task.glob("pred_*.npz"))}
    for m, f in preds.items():
        p = load_pred(f)
        assert p["queries"] == q, f"{task}/{m}: queries differ from task.tsv"
        variants = {m: p["call"]}
        if m.startswith("picrust2_") and m.removeprefix("picrust2_") in thresholds:
            variants[m + "_tuned"] = p["score"] >= thresholds[m.removeprefix("picrust2_")]
        for name, call in variants.items():
            s = ko_scores(Y, call)
            mods = scorer.matrix(pd.DataFrame(call, index=q, columns=kos)).reindex(
                columns=M_true.columns, fill_value=0.0)
            mae = (mods.to_numpy() - M_true.to_numpy()).__abs__().mean(1)
            df = pd.DataFrame(s, index=q)
            df["module_mae"] = mae
            if "copies" in p and name == m and arm == "full":
                truth_c = truth_Y.attrs["copies"].loc[q].to_numpy()
                union = (truth_c > 0) | (p["copies"] >= 0.5)
                err = np.abs(np.minimum(p["copies"], 10) - truth_c) * union
                df["copy_mae"] = err.sum(1) / np.maximum(union.sum(1), 1)
            df["method"] = name
            rows.append(df)
    out = pd.concat(rows).rename_axis("accession").reset_index()
    out[["level", "fold", "arm", "phase"]] = level, fold, arm, phase
    return out


def tune_thresholds(tasks, truth_Y):
    """PICRUSt2 presence threshold per HSP method, chosen on dev genomes (full length, all levels) to
    maximise mean per-genome F1: the same courtesy our stacker gets from dev data."""
    best = {}
    for m in P2_METHODS:
        tot, n = np.zeros(len(THRESHOLD_GRID)), 0
        for t in tasks:
            f = t / f"pred_picrust2_{m}.npz"
            info = pd.read_csv(t / "task.tsv", sep="\t")
            if not f.exists() or info.empty or info.arm.iloc[0] != "full":
                continue
            p = load_pred(f)
            Y = truth_Y.loc[p["queries"]].to_numpy()
            for i, th in enumerate(THRESHOLD_GRID):
                tot[i] += ko_scores(Y, p["score"] >= th)["ko_f1"].sum()
            n += len(p["queries"])
        if n:
            best[m] = THRESHOLD_GRID[int(np.argmax(tot))]
    return best


def pooled_scores(tasks, level, arm, stems):
    """Scores of each method for every held-out genome of one level x arm (folds stacked): (queries, {stem: S})."""
    queries, scores = [], {m: [] for m in stems}
    for t in tasks:
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if info.empty or info.level.iloc[0] != level or info.arm.iloc[0] != arm:
            continue
        if not all((t / f"pred_{m}.npz").exists() for m in stems):
            continue
        queries += info.accession.tolist()
        for m in stems:
            scores[m].append(load_pred(t / f"pred_{m}.npz")["score"])
    return queries, {m: np.vstack(v) for m, v in scores.items() if v}


def auroc_level(tasks, truth_Y, genomes, level, stems, rng, boot, min_class=5):
    """Per-KO AUROC of every method at one level (full length), PICRUSt2's best setting by mean AUROC, and the
    co-primary statistic: mean over KOs of s16fun - PICRUSt2 best, with a 95 % CI from resampling genera
    (each resample re-weights genomes by their genus's draw count and recomputes every per-KO AUROC)."""
    queries, S = pooled_scores(tasks, level, "full", stems)
    if not queries or "s16fun" not in S:
        return None, None
    Y = truth_Y.loc[queries].to_numpy()
    npos = Y.sum(0)
    keep = (npos >= min_class) & (len(Y) - npos >= min_class)
    Y = Y[:, keep]
    auc = {m: WeightedAUROC(S[m][:, keep], Y) for m in S}
    ones = np.ones(len(queries))
    A = pd.DataFrame({m: f(ones) for m, f in auc.items()}, index=truth_Y.columns[keep])
    p2 = [m for m in A if m.startswith("picrust2_")]
    best = max(p2, key=lambda m: A[m].mean())
    A["picrust2_best"] = A[best]
    codes, uniq = pd.factorize(genomes.loc[queries, "genus"])
    diffs = []
    for _ in range(boot):
        w = np.bincount(rng.integers(0, len(uniq), len(uniq)), minlength=len(uniq))[codes].astype(float)
        diffs.append(np.nanmean(auc["s16fun"](w) - auc[best](w)))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    d = (A["s16fun"] - A["picrust2_best"]).mean()
    stat = dict(auroc_s16fun=A["s16fun"].mean(), auroc_picrust2_best=A["picrust2_best"].mean(),
                auroc_picrust2_setting=best, auroc_diff=d, auroc_ci_lo=lo, auroc_ci_hi=hi, kos=int(keep.sum()),
                kos_s16fun_better=float((A["s16fun"] > A["picrust2_best"]).mean()))
    return A, stat


# ---------- analyses ----------

def add_picrust2_best(pg):
    """PICRUSt2 (best setting): per level x arm, the PICRUSt2 variant with the highest mean KO F1 on
    these genomes, chosen after seeing them (generous to PICRUSt2)."""
    p2 = pg[pg.method.str.startswith("picrust2_")]
    means = p2.groupby(["level", "arm", "method"]).ko_f1.mean()
    best = means.groupby(["level", "arm"]).idxmax().map(lambda x: x[2])
    rows = [p2[(p2.level == lv) & (p2.arm == arm) & (p2.method == m)].assign(method="picrust2_best", setting=m)
            for (lv, arm), m in best.items()]
    return pd.concat([pg] + rows, ignore_index=True), best.rename("picrust2_best_setting").reset_index()


def primary(pg, genomes, rng, boot, auroc_stats):
    """Co-primary rule (BENCHMARK.md), per level genus..class, full-length 16S, s16fun vs PICRUSt2's best setting:
    (a) per-KO AUROC: mean over KOs of the difference >= +0.05, CI lower bound > 0 (genera resampled);
    (b) KO F1: per-genome difference, CI lower bound >= -0.01 (non-inferiority, genera resampled).
    Win = (a) and (b) at every level."""
    rows = []
    for lv in PRIMARY_LEVELS:
        a = pg[(pg.level == lv) & (pg.arm == "full") & (pg.method == "s16fun")].set_index("accession").ko_f1
        b = pg[(pg.level == lv) & (pg.arm == "full") & (pg.method == "picrust2_best")].set_index("accession").ko_f1
        common = a.index.intersection(b.index)
        if common.empty or lv not in auroc_stats:
            continue
        d = a[common] - b[common]
        mean, lo, hi = cluster_bootstrap(d.to_numpy(), genomes.loc[common, "genus"].to_numpy(), rng, boot)
        st = auroc_stats[lv]
        auroc_ok = bool(st["auroc_diff"] >= AUROC_MARGIN and st["auroc_ci_lo"] > 0)
        f1_ok = bool(lo >= -F1_NONINFERIORITY)
        rows.append(dict(level=lv, genomes=len(common), genera=genomes.loc[common, "genus"].nunique(),
                         **{k: v for k, v in st.items()}, auroc_superior=auroc_ok,
                         f1_s16fun=a[common].mean(), f1_picrust2_best=b[common].mean(), f1_diff=mean, f1_ci_lo=lo,
                         f1_ci_hi=hi, f1_noninferior=f1_ok, genomes_s16fun_better=(d > 0).mean(),
                         genomes_picrust2_better=(d < 0).mean(), win=auroc_ok and f1_ok))
    return pd.DataFrame(rows)


def secondary(pg, genomes, rng, boot):
    """Paired s16fun - PICRUSt2 best for the other endpoints (positive = s16fun better), Holm-adjusted
    over the endpoint x level family by the bootstrap's two-sided p (share of resamples crossing 0)."""
    rows = []
    sign = {"module_mae": -1, "copy_mae": -1}
    for metric in ["ko_precision", "ko_recall", "module_mae", "copy_mae"]:
        for lv in LEVELS:
            for arm in sorted(pg.arm.unique()):
                a = pg[(pg.level == lv) & (pg.arm == arm) & (pg.method == "s16fun")].set_index("accession")
                b = pg[(pg.level == lv) & (pg.arm == arm) & (pg.method == "picrust2_best")].set_index("accession")
                if metric not in a or metric not in b:
                    continue
                common = a.index.intersection(b.index)
                d = (a.loc[common, metric] - b.loc[common, metric]).dropna() * sign.get(metric, 1)
                if d.empty:
                    continue
                g = genomes.loc[d.index, "genus"].to_numpy()
                mean, lo, hi = cluster_bootstrap(d.to_numpy(), g, rng, boot)
                rows.append(dict(endpoint=metric, level=lv, arm=arm, genomes=len(d), s16fun=a.loc[d.index, metric].mean(),
                                 picrust2_best=b.loc[d.index, metric].mean(), diff_s16fun_better=mean, ci_lo=lo, ci_hi=hi))
    return pd.DataFrame(rows)


def novelty_table(pg, nov):
    d = pg[(pg.arm == "full") & pg.method.isin(["s16fun", "picrust2_best"])].merge(
        nov, on=["accession", "level", "fold"], how="left")
    d["identity_bin"] = pd.cut(d.top_id, [0, 0.85, 0.90, 0.94, 0.97, 0.99, 1.0001],
                               labels=["<85 %", "85-90 %", "90-94 %", "94-97 %", "97-99 %", ">=99 %"])
    d["nsti_bin"] = pd.cut(d.nsti, [-1, 0.03, 0.08, 0.16, 0.36, 100],
                           labels=["<0.03", "0.03-0.08", "0.08-0.16", "0.16-0.36", ">0.36"])
    rows = []
    for by in ["identity_bin", "nsti_bin"]:
        t = d.groupby([by, "method"], observed=True).agg(ko_f1=("ko_f1", "mean"), module_mae=("module_mae", "mean"),
                                                          n=("ko_f1", "size")).reset_index().rename(columns={by: "bin"})
        t["bin_by"] = by
        rows.append(t)
    return pd.concat(rows), d


def resources(traces):
    if not traces:
        return pd.DataFrame()
    t = pd.concat(pd.read_csv(f, sep="\t") for f in traces)
    t = t[t.status.isin(["COMPLETED", "CACHED"])].drop_duplicates("name", keep="last")
    t["process"] = t.name.str.replace(r" \(.*\)$", "", regex=True)

    def to_s(x):
        unit = {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}
        try:
            return sum(float(p[:-2]) * unit["ms"] if p.endswith("ms") else float(p[:-1]) * unit[p[-1]]
                       for p in str(x).split())
        except (ValueError, KeyError):  # "-" or an unexpected format
            return np.nan

    def to_gb(x):
        parts = str(x).split()
        if len(parts) != 2:            # "-", "nan" or a bare "0" (no memory recorded)
            return float(parts[0]) * 1e-9 if parts and parts[0].replace(".", "", 1).isdigit() else np.nan
        return float(parts[0]) * {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1, "TB": 1e3}[parts[1]]
    t["realtime_s"] = t.realtime.map(to_s)
    t["cpu_h"] = t.realtime_s * t["%cpu"].astype(str).str.rstrip("%").astype(float) / 100 / 3600
    t["peak_gb"] = t.peak_rss.map(to_gb)
    return t.groupby("process").agg(tasks=("name", "size"), wall_h_sum=("realtime_s", lambda s: s.sum() / 3600),
                                    cpu_h=("cpu_h", "sum"), median_task_min=("realtime_s", lambda s: s.median() / 60),
                                    peak_gb=("peak_gb", "max")).reset_index()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--tasks", nargs="+", required=True)
    ap.add_argument("--phase", required=True, choices=["dev", "confirm"])
    ap.add_argument("--thresholds")
    ap.add_argument("--stacker-importance")
    ap.add_argument("--community", help="bench_community.py output dir")
    ap.add_argument("--trace", nargs="*", default=[])
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--auroc-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpus", type=int, default=16)
    ap.add_argument("--max-ordination", type=int, default=1500)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    genomes, kos, copies, _ = load_data(a.data)
    truth_Y = pd.DataFrame((copies > 0).toarray(), index=genomes.index, columns=kos)
    truth_Y.attrs["copies"] = pd.DataFrame(np.minimum(copies.toarray(), 10).astype(np.float32),
                                           index=genomes.index, columns=kos)
    scorer = ModuleScorer()
    truth_mod = scorer.matrix(truth_Y)
    tasks = sorted(Path(t) for t in a.tasks)

    # leak and consistency checks across tasks
    splits = pd.read_csv(a.splits, sep="\t")
    for t in tasks:
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if info.empty:
            continue
        lv, fold = info.level.iloc[0], info.fold.iloc[0]
        held = set(splits[(splits.level == lv) & (splits.fold == fold)].accession)
        train = set((t / "train.txt").read_text().split())
        assert set(info.accession) <= held and not set(info.accession) & train, f"{t}: query/training overlap"
        assert (info.phase == a.phase).all(), f"{t}: phase mismatch"
        part = splits.drop_duplicates("accession").set_index("accession").partition
        assert (part.loc[info.accession] == a.phase).all(), f"{t}: genome from the other partition"

    if a.phase == "dev":
        thresholds = tune_thresholds(tasks, truth_Y)
        pd.Series(thresholds, name="threshold").rename_axis("hsp_method").to_csv(out / "picrust2_thresholds.tsv", sep="\t")
    else:
        thresholds = pd.read_csv(a.thresholds, sep="\t", index_col=0)["threshold"].to_dict() if a.thresholds else {}

    pg = pd.concat(Parallel(n_jobs=a.cpus)(delayed(task_metrics)(t, truth_Y, truth_mod, kos, thresholds, scorer)
                                           for t in tasks), ignore_index=True)
    pg, best_setting = add_picrust2_best(pg)
    pg.to_parquet(out / "per_genome.parquet")
    best_setting.to_csv(out / "picrust2_best_setting.tsv", sep="\t", index=False)

    nov = pd.concat(pd.read_csv(t / "novelty.tsv", sep="\t").assign(**pd.read_csv(t / "task.tsv", sep="\t").iloc[0][
        ["level", "fold", "arm"]].to_dict()) for t in tasks if (t / "novelty.tsv").exists()).rename(columns={"query": "accession"})
    nsti = pd.concat(pd.read_csv(t / "nsti.tsv", sep="\t").assign(**pd.read_csv(t / "task.tsv", sep="\t").iloc[0][
        ["level", "fold", "arm"]].to_dict()) for t in tasks if (t / "nsti.tsv").exists()).rename(columns={"query": "accession"})
    nov = nov[nov.arm == "full"].merge(nsti[nsti.arm == "full"][["accession", "level", "fold", "nsti"]],
                                       on=["accession", "level", "fold"], how="left").drop(columns="arm")

    summary = pg.groupby(["level", "arm", "method"])[["ko_f1", "ko_precision", "ko_recall", "module_mae"]].mean().reset_index()
    summary.to_csv(out / "summary.tsv", sep="\t", index=False)
    stems = ["s16fun"] + [f"picrust2_{m}" for m in P2_METHODS]
    res_auc = Parallel(n_jobs=min(len(PRIMARY_LEVELS), a.cpus))(
        delayed(auroc_level)(tasks, truth_Y, genomes, lv, stems, np.random.default_rng(a.seed + i), a.auroc_boot)
        for i, lv in enumerate(PRIMARY_LEVELS))
    auroc = {lv: A for lv, (A, _) in zip(PRIMARY_LEVELS, res_auc) if A is not None}
    auroc_stats = {lv: st for lv, (_, st) in zip(PRIMARY_LEVELS, res_auc) if st is not None}
    if auroc:
        pd.concat(auroc, names=["level", "ko"]).to_csv(out / "per_ko_auroc.tsv.gz", sep="\t")
    prim = primary(pg, genomes, rng, a.boot, auroc_stats)
    prim.to_csv(out / "primary.tsv", sep="\t", index=False)
    sec = secondary(pg, genomes, rng, a.boot)
    sec.to_csv(out / "secondary.tsv", sep="\t", index=False)
    novt, novd = novelty_table(pg, nov)
    novt.to_csv(out / "by_novelty.tsv", sep="\t", index=False)
    res = resources(a.trace)
    res.to_csv(out / "resources.tsv", sep="\t", index=False)

    # figures
    fig_methods_diagram(out / "fig0_methods.png")
    fig_primary(summary, prim, out / "fig1_primary.png")
    fig_novelty(novt, out / "fig2_novelty.png")
    fig_phyla(pg, genomes, out / "fig3_phyla.png")
    if "genus" in auroc:
        fig_per_ko(auroc["genus"], out / "fig4_per_ko.png")
    fig_arms(summary, out / "fig5_arms.png")
    ordi = None
    gb = best_setting[(best_setting.level == "genus") & (best_setting.arm == "full")].picrust2_best_setting
    if len(gb):
        stem = gb.iloc[0].removeprefix("pred_").removesuffix("_tuned")
        th = thresholds.get(stem.removeprefix("picrust2_")) if gb.iloc[0].endswith("_tuned") else None
        ordi = ordination(tasks, truth_Y, truth_mod, scorer, rng, a.max_ordination, "genus",
                          {"s16fun": ("s16fun", None), "picrust2_best": (stem, th), "knn5": ("knn5", None),
                           "kmer_mlp": ("kmer_mlp", None)})
    if ordi is not None:
        fig_ordination(ordi, out / "fig6_ordination.png")
        ordi["stats"].to_csv(out / "ordination_stats.tsv", sep="\t", index=False)
    if not res.empty:
        fig_resources(res, out / "fig7_resources.png")

    community = None
    if a.community and (Path(a.community) / "community_summary.tsv").exists():
        import shutil
        for f in ["fig_community.png", "community_summary.tsv", "community_reps.tsv"]:
            shutil.copy(Path(a.community) / f, out / ("fig8_community.png" if f.endswith(".png") else f))
        community = pd.read_csv(out / "community_summary.tsv", sep="\t")

    from bench_report import write_report
    write_report(out, a.phase, genomes, summary, prim, sec, novt, best_setting, thresholds, res, ordi,
                 a.stacker_importance, pg, community)


if __name__ == "__main__":
    main()
