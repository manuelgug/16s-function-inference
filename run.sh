#!/usr/bin/env bash
# One command runs the whole project: downloads tools and databases, then the Nextflow pipeline.
#
#   ./run.sh bench_test                   # benchmark smoke test: ~1,500 genomes of PICRUSt2's reference
#   ./run.sh bench                        # benchmark: s16fun vs PICRUSt2 and baselines (BENCH_PHASE dev|confirm)
#   ./run.sh test                         # GTDB genome pipeline (annotation, Aim 2): 24 genomes
#   ./run.sh pilot                        # ~2,000 genomes
#   ./run.sh full                         # every genome passing the filters
#   ./run.sh pilot --min_ssu_len 1300     # any parameter below can be overridden per run
#   CONTAINER=docker ./run.sh test        # use Docker instead of Apptainer
#
# Interrupted? Run the same command again: downloads resume and Nextflow skips finished work.
set -euo pipefail
cd "$(dirname "$0")"

# ============================== USER PARAMETERS ==============================
# Edit here. Command-line flags (--name value) override these for a single run.
# Changing a value re-runs only the steps that depend on it.

# --- Reference data versions ---
GTDB_RELEASE=232            # GTDB release (DECISIONS D013)
SILVA_RELEASE=144           # SILVA SSU Ref NR99 release for SINA (D023)

# --- Genome filters (D014) ---
MIN_COMPLETENESS=90         # CheckM2 completeness, %
MAX_CONTAMINATION=5         # CheckM2 contamination, %
MIN_SSU_LEN=1400            # minimum 16S length, bp (GTDB prefilter and our barrnap copies)

# --- Sampling (full scale only). Pilot and test use the fixed lists committed in config/, so the
#     genome filters above change them only if you regenerate the lists (see bin/select_genomes.py). ---
FULL_N=0                    # 0 = every genome passing the filters; >0 = stratified sample of that size
SEED=42
MIN_PER_PHYLUM=10

# --- 16S features (D017) ---
MIN_COL_COVERAGE=0.95       # keep alignment columns with a base in >= this fraction of genomes
MIN_ALLELE_FREQ=0.01        # keep allele features with frequency in [f, 1-f]

# --- In-silico PCR primers (D015): name, forward 5'->3', reverse 5'->3'. IUPAC codes allowed.
# Add, remove or edit lines freely; each line becomes one amplicon arm.
PRIMERS="
V4    GTGYCAGCMGCCGCGGTAA   GGACTACNVGGGTWTCTAAT     515F-Y / 806RB
V3V4  CCTACGGGNGGCWGCAG     GACTACHVGGGTATCTAATCC    341F / 805R
"
PCR_ERROR_RATE=0.1          # cutadapt max error rate per primer

# --- Function targets (D017), used by Aim 2 ---
KO_MIN_PREVALENCE=0.05
KO_MAX_PREVALENCE=0.95

# --- Benchmark (DESIGN §B, BENCHMARK.md) ---
HOLDOUT_LEVELS="random,genus,family,order,class"   # whole taxa held out per fold (random = species)
CV_FOLDS=5
BENCH_FOLDS="0,1,2,3,4"     # folds to run (e.g. "0" for exploration; the rest reuse its cached tasks later)
BENCH_PHASE=dev             # dev: develop + tune on dev families | confirm: frozen settings, confirm families, once
BENCH_ARMS="full,V4,V3V4"   # full-length 16S and the amplicons named in PRIMERS
BENCH_P2_METHODS="mp,emp_prob,pic,subtree_average"   # PICRUSt2 hidden-state prediction settings, all run
BENCH_MIN_KO_GENOMES=10     # KOs scored: present in >= this many reference genomes
BENCH_DEV_FRACTION=0.5      # share of families in the dev partition (the rest: confirm)
BENCH_SUBSAMPLE=1500        # bench_test only: genomes kept (bench uses all 27,870)
BENCH_STACKER_GENOMES=150   # held-out genomes per dev task sampled to train the stacker
BENCH_BOOT=2000             # cluster-bootstrap resamples (genera)
BENCH_FROZEN="$PWD/config/bench_frozen"   # stackers + PICRUSt2 thresholds frozen after dev (confirm reads them)

# --- Compute ---
BATCH_SIZE=100              # genomes per per-genome task
MAX_CPUS=$(nproc)
MAX_MEMORY="$(awk '/MemTotal/{printf "%d", $2/1048576*0.9}' /proc/meminfo).GB"   # 90% of RAM
DB_DIR="${DB_DIR:-$PWD/db}"  # ~65 GB of databases; set DB_DIR in the environment to share one copy
OUTDIR="$PWD/results"
# =============================================================================

NXF_VERSION=26.10.4
JAVA_MIN=17
CONTAINER=${CONTAINER:-apptainer}

die() { echo "ERROR: $*" >&2; exit 1; }
PROFILE=${1:-}
[[ "$PROFILE" =~ ^(test|pilot|full|bench_test|bench)$ ]] || die "usage: ./run.sh {bench_test|bench|test|pilot|full} [--param value ...]"
shift
[[ "$CONTAINER" =~ ^(apptainer|docker)$ ]] || die "CONTAINER must be apptainer or docker"
command -v "$CONTAINER" >/dev/null || die "$CONTAINER not found. Install Apptainer (https://apptainer.org/docs/admin/main/installation.html) or Docker, or set CONTAINER=docker"
for t in curl gzip tar md5sum; do command -v $t >/dev/null || die "$t not found"; done
mkdir -p .bin

# --- Java (Nextflow needs >= 17): use the system one, or fetch a private Temurin JRE into .bin/ ---
java_ok() { "$1" -version 2>&1 | awk -F'"' '/version/{split($2,v,"."); exit !(v[1] >= '"$JAVA_MIN"')}'; }
if command -v java >/dev/null && java_ok java; then
  JAVA=java
else
  if [ ! -x .bin/jre/bin/java ]; then
    echo "[run.sh] installing Temurin JRE 21 into .bin/jre"
    curl -fL --retry 5 -o .bin/jre.tar.gz "https://api.adoptium.net/v3/binary/latest/21/ga/linux/$(uname -m | sed 's/x86_64/x64/')/jre/hotspot/normal/eclipse"
    mkdir -p .bin/jre && tar -xzf .bin/jre.tar.gz -C .bin/jre --strip-components=1 && rm .bin/jre.tar.gz
  fi
  JAVA="$PWD/.bin/jre/bin/java"
fi
export JAVA_CMD="$JAVA"

# --- Nextflow (pinned, self-contained jar) ---
NXF=".bin/nextflow-$NXF_VERSION"
if [ ! -s "$NXF" ]; then
  echo "[run.sh] installing Nextflow $NXF_VERSION into .bin/"
  curl -fL --retry 5 -o "$NXF.part" "https://github.com/nextflow-io/nextflow/releases/download/v$NXF_VERSION/nextflow-$NXF_VERSION-dist"
  chmod +x "$NXF.part" && mv "$NXF.part" "$NXF"
fi

# --- Databases (the benchmark needs none: PICRUSt2's reference ships in its container) ---
[[ "$PROFILE" == bench* ]] || scripts/setup_db.sh "$DB_DIR" "$GTDB_RELEASE" "$SILVA_RELEASE"
[[ "$PROFILE" == bench_test ]] || BENCH_SUBSAMPLE=0

# --- Parameters -> files Nextflow reads ---
RUN_DIR="$OUTDIR/$PROFILE"
mkdir -p "$RUN_DIR"
printf '%s\n' "$PRIMERS" | awk 'NF >= 3 {print $1 "\t" $2 "\t" $3}' > "$RUN_DIR/primers.tsv"
[ -s "$RUN_DIR/primers.tsv" ] || die "PRIMERS is empty"
cat > "$RUN_DIR/params.json" <<EOF
{
  "gtdb_release": $GTDB_RELEASE,
  "silva_release": "$SILVA_RELEASE",
  "min_completeness": $MIN_COMPLETENESS,
  "max_contamination": $MAX_CONTAMINATION,
  "min_ssu_len": $MIN_SSU_LEN,
  "full_n": $FULL_N,
  "seed": $SEED,
  "min_per_phylum": $MIN_PER_PHYLUM,
  "min_col_coverage": $MIN_COL_COVERAGE,
  "min_allele_freq": $MIN_ALLELE_FREQ,
  "primers": "$RUN_DIR/primers.tsv",
  "pcr_error_rate": $PCR_ERROR_RATE,
  "ko_min_prevalence": $KO_MIN_PREVALENCE,
  "ko_max_prevalence": $KO_MAX_PREVALENCE,
  "holdout_levels": "$HOLDOUT_LEVELS",
  "cv_folds": $CV_FOLDS,
  "bench_phase": "$BENCH_PHASE",
  "bench_arms": "$BENCH_ARMS",
  "bench_p2_methods": "$BENCH_P2_METHODS",
  "bench_min_ko_genomes": $BENCH_MIN_KO_GENOMES,
  "bench_dev_fraction": $BENCH_DEV_FRACTION,
  "bench_subsample": $BENCH_SUBSAMPLE,
  "bench_stacker_genomes": $BENCH_STACKER_GENOMES,
  "bench_boot": $BENCH_BOOT,
  "bench_frozen": "$BENCH_FROZEN",
  "bench_folds": "$BENCH_FOLDS",
  "batch_size": $BATCH_SIZE,
  "max_cpus": $MAX_CPUS,
  "max_memory": "$MAX_MEMORY",
  "db_dir": "$DB_DIR",
  "outdir": "$RUN_DIR"
}
EOF

echo "[run.sh] profile=$PROFILE container=$CONTAINER outdir=$RUN_DIR"
# Launch from the profile's own directory so -resume always resumes that profile's last run.
PROJECT=$PWD
cd "$RUN_DIR"
SCRIPT=main.nf
[[ "$PROFILE" == bench* ]] && SCRIPT=benchmark.nf
exec "$PROJECT/$NXF" run "$PROJECT/$SCRIPT" -profile "$PROFILE,$CONTAINER" -params-file params.json \
  -work-dir "$PROJECT/work/$PROFILE" -resume "$@"
