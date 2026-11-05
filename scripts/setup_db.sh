#!/usr/bin/env bash
# Download, verify and unpack the reference databases into DB_DIR.
# Idempotent: finished files carry a .ok marker and are skipped.
# Resumable: interrupted downloads continue where they stopped (curl -C -).
# Host needs only: curl, gzip, tar, md5sum.
#
# Usage: scripts/setup_db.sh DB_DIR GTDB_RELEASE SILVA_RELEASE
set -euo pipefail

DB_DIR=${1:?usage: setup_db.sh DB_DIR GTDB_RELEASE SILVA_RELEASE}
GTDB=${2:?GTDB release, e.g. 232}
SILVA=${3:?SILVA release, e.g. 144}
HERE=$(cd "$(dirname "$0")/.." && pwd)
CHECKSUMS="$HERE/config/db_checksums.tsv"   # md5s for sources that publish none (trust on first use)
mkdir -p "$DB_DIR"
MANIFEST="$DB_DIR/MANIFEST.tsv"
[ -s "$MANIFEST" ] || printf 'file\tsource_url\tmd5\tdate\n' > "$MANIFEST"

log() { printf '[setup_db] %s\n' "$*" >&2; }

# fetch URL DEST [EXPECTED_MD5]: download to DEST, verify, mark DEST.ok
fetch() {
  local url=$1 dest=$2 want=${3:-}
  [ -e "$dest.ok" ] && return 0
  mkdir -p "$(dirname "$dest")"
  log "downloading $url"
  curl -fL --retry 10 --retry-delay 30 -C - -o "$dest.part" "$url"
  local got; got=$(md5sum "$dest.part" | cut -d' ' -f1)
  if [ -z "$want" ]; then
    want=$(awk -v u="$url" -F'\t' '$1==u{print $2}' "$CHECKSUMS" 2>/dev/null || true)
  fi
  if [ -n "$want" ] && [ "$got" != "$want" ]; then
    log "MD5 MISMATCH for $url: got $got, expected $want. Deleting partial file; rerun to retry."
    rm -f "$dest.part"; exit 1
  fi
  [ -n "$want" ] || log "WARNING: no published md5 for $url; add this line to config/db_checksums.tsv: $url	$got"
  mv "$dest.part" "$dest"
  printf '%s\t%s\t%s\t%s\n' "${dest#"$DB_DIR"/}" "$url" "$got" "$(date -u +%FT%TZ)" >> "$MANIFEST"
  touch "$dest.ok"
}

# unpack ARCHIVE OUTDIR MARKER: gunzip or untar once, then delete the archive (its .ok marker stays)
unpack() {
  local arc=$1 out=$2 marker=$3
  [ -e "$marker" ] && return 0
  log "unpacking $arc"
  case "$arc" in
    *.tar.gz) tar -xzf "$arc" -C "$out" ;;
    *.gz)     gunzip -c "$arc" > "${arc%.gz}.part" && mv "${arc%.gz}.part" "${arc%.gz}" ;;
  esac
  touch "$marker"
  rm -f "$arc"
}

# --- GTDB: metadata, taxonomy, trees (md5s from the release's MD5SUM.txt) ---
G="$DB_DIR/gtdb_r$GTDB"; mkdir -p "$G"
GURL="https://data.gtdb.ecogenomic.org/releases/release$GTDB/$GTDB.0"
GMD5=$(curl -fsL --retry 5 "$GURL/MD5SUM.txt")
for f in bac120_metadata_r$GTDB.tsv.gz ar53_metadata_r$GTDB.tsv.gz \
         bac120_taxonomy_r$GTDB.tsv ar53_taxonomy_r$GTDB.tsv \
         bac120_r$GTDB.tree ar53_r$GTDB.tree; do
  fetch "$GURL/$f" "$G/$f" "$(awk -v f="./$f" '$2==f{print $1}' <<<"$GMD5")"
done

# --- SILVA SSU Ref NR99 ARB file for SINA (md5 published next to it) ---
S="$DB_DIR/silva_$SILVA"; mkdir -p "$S"
SF="SILVA_${SILVA}_SSURef_NR99_opt.arb.gz"
SURL="https://ftp.arb-silva.de/release_$SILVA/ARB_files/$SF"
if [ ! -e "$S/unpacked.ok" ]; then
  fetch "$SURL" "$S/$SF" "$(curl -fsL --retry 5 "$SURL.md5" | cut -d' ' -f1)"
  unpack "$S/$SF" "$S" "$S/unpacked.ok"
fi

# --- eggNOG-mapper data, eggNOG 5.0.2 (no published md5: checked against config/db_checksums.tsv) ---
E="$DB_DIR/eggnog_5.0.2"; mkdir -p "$E"
EURL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
for f in eggnog.db.gz eggnog_proteins.dmnd.gz eggnog.taxa.tar.gz; do
  if [ ! -e "$E/$f.unpacked.ok" ]; then
    fetch "$EURL/$f" "$E/$f"
    unpack "$E/$f" "$E" "$E/$f.unpacked.ok"
  fi
done

touch "$DB_DIR/setup.ok"
log "all databases ready in $DB_DIR"
