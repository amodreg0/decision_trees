#!/bin/bash
#SBATCH --job-name=truvari
#SBATCH --output=truvari_%A_%a.out
#SBATCH --error=truvari_%A_%a.err
#SBATCH --array=0-7
#SBATCH --time=00:30:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5
#SBATCH --nodelist=node

set -euo pipefail
: "${PS1:=}"

eval "$(conda shell.bash hook)"
conda activate truvari

# ----------------------------
# Filtering thresholds from decision tree:
#
# ALL — pruned_min:
#
# Keep if:
#
#   1) SOMATICSCORE >= 50 AND SVTYPE in {DEL,DUP,INV}
#
#   OR
#
#   2) SOMATICSCORE >= 59 AND SVTYPE in {BND,INS}
#      AND ND >= 2.9
#      AND PR_SR_ratio = PR_alt / (PR_alt + SR_alt) >= 0.44
#
# Remove otherwise.
#
# Notes:
#   - SS in the tree corresponds to INFO/SOMATICSCORE.
#   - ND is read from INFO/ND2 first, then INFO/ND if present.
#   - PR/SR are FORMAT fields from the tumor sample.
#   - Missing/unparseable SOMATICSCORE removes the variant.
#   - Missing/unparseable SVTYPE removes the variant.
#   - For BND/INS, missing/unparseable ND removes the variant.
#   - For BND/INS, missing/unparseable PR/SR removes the variant.
# ----------------------------
MIN_SOMATICSCORE_MAIN=50
MIN_SOMATICSCORE_BND_INS=59
MIN_ND=2.9
MIN_PR_SR_RATIO=0.44

# Tumor sample index for PR/SR fields, 0-based among samples in VCF
TUMOR_SAMPLE_INDEX=1

# ----------------------------
# Paths
# ----------------------------
SEVERUS_BASE_PATH="../SVs"
SEVERUS_GROUND_TRUTH_BASE_PATH="../SVs/severus"
TRUVARI_OUTPUT_BASE="../SVs/truvari_tree_filters"
REFERENCE_PATH="../hs1.fa"

CELL_LINES=("H1437" "H2009" "HCC1954" "HCC1395" "HCC1937" "Hs578T" "HG008" "COLO829")
CELL_LINE="${CELL_LINES[$SLURM_ARRAY_TASK_ID]}"

echo "Processing cell line: $CELL_LINE"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"

GROUND_TRUTH="$SEVERUS_GROUND_TRUTH_BASE_PATH/${CELL_LINE}_merged/tp-base.clean.vcf.gz"

# ND + in_bed annotated VCFs
SEVERUS_HG38="$SEVERUS_BASE_PATH/$CELL_LINE/Severus_${CELL_LINE}_hg38_lifted_ND_annot.vcf.gz"
SEVERUS_T2T="$SEVERUS_BASE_PATH/$CELL_LINE/Severus_${CELL_LINE}_t2t_ND_annot.vcf.gz"

mkdir -p "$TRUVARI_OUTPUT_BASE/${CELL_LINE}"

HG38_T2T="$TRUVARI_OUTPUT_BASE/${CELL_LINE}/HG38_T2T"
HG38_GROUNDTRUTH="$TRUVARI_OUTPUT_BASE/${CELL_LINE}/HG38_GT"
T2T_GROUNDTRUTH="$TRUVARI_OUTPUT_BASE/${CELL_LINE}/T2T_GT"
ALL3="$TRUVARI_OUTPUT_BASE/${CELL_LINE}/ALL3"

# ----------------------------
# Temp dir
# ----------------------------
TMPROOT="${SLURM_TMPDIR:-$TRUVARI_OUTPUT_BASE/${CELL_LINE}}"
TMPDIR="$(mktemp -d -p "$TMPROOT" "truvari_filter_${CELL_LINE}_XXXXXX")"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

# ----------------------------
# Requirements
# ----------------------------
for cmd in python3 bgzip tabix truvari; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "ERROR: missing required command: $cmd" >&2
    exit 1
  }
done

# ----------------------------
# Filter Severus VCF using decision-tree rules
# ----------------------------
filter_severus_streaming() {
  local in_vcf="$1"
  local out_vcf="$2"

  python3 - \
    "$in_vcf" \
    "${MIN_SOMATICSCORE_MAIN}" \
    "${MIN_SOMATICSCORE_BND_INS}" \
    "${MIN_ND}" \
    "${MIN_PR_SR_RATIO}" \
    "${TUMOR_SAMPLE_INDEX}" \
  <<'PY' \
  | bgzip -c > "$out_vcf"
import sys, gzip

in_vcf                = sys.argv[1]
min_somatic_main      = float(sys.argv[2])
min_somatic_bnd_ins   = float(sys.argv[3])
min_nd                = float(sys.argv[4])
min_pr_sr_ratio       = float(sys.argv[5])
tumor_idx             = int(sys.argv[6])

def parse_info(info_str):
    d = {}
    if info_str == ".":
        return d
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" not in item:
            # Flag-style INFO field
            d[item] = True
            continue
        k, v = item.split("=", 1)
        d[k] = v
    return d

def parse_format(fmt, sample):
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: (vals[i] if i < len(vals) else ".") for i, k in enumerate(keys)}

def parse_float(v):
    if v is None or v in (".", ""):
        return None
    try:
        return float(v)
    except ValueError:
        return None

def parse_allele_count(v):
    """
    Parse a comma-separated ref,alt FORMAT field, e.g. PR or SR.
    Returns (alt_count, ok).
    ok=False when field is absent, missing, or unparseable.
    """
    if v is None or v in (".", ""):
        return (0, False)
    parts = v.split(",")
    if len(parts) < 2:
        return (0, False)
    try:
        return (int(parts[1]), True)
    except ValueError:
        return (0, False)

def normalize_svtype(svtype):
    if svtype is None or svtype in (".", ""):
        return None

    # Sometimes SVTYPE can theoretically contain commas. Use the first value.
    svtype = svtype.split(",")[0].strip().upper()

    # Normalize common aliases if needed.
    if svtype == "DUP:TANDEM":
        svtype = "DUP"

    return svtype

sample_count = None
header_has_pr = False
header_has_sr = False
header_has_som = False
header_has_svtype = False
header_has_nd = False
header_has_nd2 = False

opener = gzip.open if in_vcf.endswith(".gz") else open

for line in opener(in_vcf, "rt"):

    # ---- meta-information lines ----
    if line.startswith("##"):
        if line.startswith("##FORMAT=<ID=PR"):
            header_has_pr = True
        elif line.startswith("##FORMAT=<ID=SR"):
            header_has_sr = True
        elif line.startswith("##INFO=<ID=SOMATICSCORE"):
            header_has_som = True
        elif line.startswith("##INFO=<ID=SVTYPE"):
            header_has_svtype = True
        elif line.startswith("##INFO=<ID=ND,"):
            header_has_nd = True
        elif line.startswith("##INFO=<ID=ND2"):
            header_has_nd2 = True

        sys.stdout.write(line)
        continue

    # ---- #CHROM header line ----
    if line.startswith("#CHROM"):
        cols = line.rstrip("\n").split("\t")
        sample_count = max(0, len(cols) - 9)

        if sample_count == 0:
            raise SystemExit("ERROR: VCF has no sample columns")

        if tumor_idx < 0 or tumor_idx >= sample_count:
            raise SystemExit(
                f"ERROR: TUMOR_SAMPLE_INDEX={tumor_idx} out of range for {sample_count} samples"
            )

        if not header_has_som:
            raise SystemExit("ERROR: INFO/SOMATICSCORE not present in VCF header")

        if not header_has_svtype:
            raise SystemExit("ERROR: INFO/SVTYPE not present in VCF header")

        if not header_has_nd and not header_has_nd2:
            raise SystemExit("ERROR: neither INFO/ND nor INFO/ND2 present in VCF header")

        if not header_has_pr or not header_has_sr:
            raise SystemExit("ERROR: FORMAT/PR and/or FORMAT/SR not present in VCF header")

        sys.stdout.write(line)
        continue

    # ---- pass through any remaining header / blank lines ----
    if not line or line.startswith("#"):
        sys.stdout.write(line)
        continue

    # ---- variant record ----
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 10:
        # Cannot evaluate the tree without FORMAT/sample fields.
        continue

    info = parse_info(cols[7])

    # ------------------------------------------------------------------
    # Node 1:
    #   SS < 50    -> FALSE/remove
    #   SS >= 50   -> continue
    #
    # SS corresponds to INFO/SOMATICSCORE.
    # Missing/unparseable SOMATICSCORE -> remove.
    # ------------------------------------------------------------------
    somatic = parse_float(info.get("SOMATICSCORE"))

    if somatic is None:
        continue

    if somatic < min_somatic_main:
        continue

    # ------------------------------------------------------------------
    # Node 2:
    #   SVTYPE = DEL,DUP,INV -> TRUE/keep
    #   SVTYPE = BND,INS     -> continue
    #   Other/missing SVTYPE -> remove
    # ------------------------------------------------------------------
    svtype = normalize_svtype(info.get("SVTYPE"))

    if svtype in {"DEL", "DUP", "INV"}:
        sys.stdout.write(line)
        continue

    if svtype not in {"BND", "INS"}:
        continue

    # ------------------------------------------------------------------
    # For BND/INS only:
    #
    # Node 3:
    #   ND < 2.9   -> FALSE/remove
    #   ND >= 2.9  -> continue
    #
    # Prefer INFO/ND2 if present, otherwise INFO/ND.
    # Missing/unparseable ND -> remove.
    # ------------------------------------------------------------------
    nd = parse_float(info.get("ND2", info.get("ND")))

    if nd is None:
        continue

    if nd < min_nd:
        continue

    # ------------------------------------------------------------------
    # Node 4:
    #   PR_SR_ratio < 0.44   -> FALSE/remove
    #   PR_SR_ratio >= 0.44  -> continue
    #
    # PR_SR_ratio = PR_alt / (PR_alt + SR_alt)
    # Missing/unparseable PR or SR -> remove.
    # ------------------------------------------------------------------
    if len(cols) < 9 + sample_count:
        continue

    fmt = cols[8]
    tumor_sample = cols[9 + tumor_idx]
    fm = parse_format(fmt, tumor_sample)

    pr_alt, pr_ok = parse_allele_count(fm.get("PR"))
    sr_alt, sr_ok = parse_allele_count(fm.get("SR"))

    if not pr_ok or not sr_ok:
        continue

    total = pr_alt + sr_alt
    if total == 0:
        continue

    pr_sr_ratio = pr_alt / float(total)

    if pr_sr_ratio < min_pr_sr_ratio:
        continue

    # ------------------------------------------------------------------
    # Node 5:
    #   SS < 59    -> FALSE/remove
    #   SS >= 59   -> TRUE/keep
    #
    # This applies only to BND/INS after ND and PR_SR_ratio pass.
    # ------------------------------------------------------------------
    if somatic < min_somatic_bnd_ins:
        continue

    sys.stdout.write(line)
PY

  tabix -f -p vcf "$out_vcf"

  local n
  n="$(python3 - <<PY
import gzip
c = 0
with gzip.open("$out_vcf", "rt") as f:
    for line in f:
        if line and not line.startswith("#"):
            c += 1
print(c)
PY
)"
  echo "Filtered according to decision tree -> $out_vcf ; variants kept: $n"
}

# ----------------------------
# GT VCF: pass through unchanged
# No somatic/ND/PR-SR filters on ground truth.
# ----------------------------
passthrough_gt() {
  local in_vcf="$1"
  local out_vcf="$2"

  cp "$in_vcf" "$out_vcf"
  tabix -f -p vcf "$out_vcf"

  local n
  n="$(python3 - <<PY
import gzip
c = 0
with gzip.open("$out_vcf", "rt") as f:
    for line in f:
        if line and not line.startswith("#"):
            c += 1
print(c)
PY
)"
  echo "Pass-through GT -> $out_vcf ; variants: $n"
}

# ----------------------------
# Build filtered inputs
# ----------------------------
FILT_HG38="$TMPDIR/Severus_${CELL_LINE}_hg38.filtered.vcf.gz"
FILT_T2T="$TMPDIR/Severus_${CELL_LINE}_t2t.filtered.vcf.gz"
FILT_GT="$TMPDIR/${CELL_LINE}_GT.filtered.vcf.gz"

filter_severus_streaming "$SEVERUS_HG38" "$FILT_HG38"
filter_severus_streaming "$SEVERUS_T2T"  "$FILT_T2T"
passthrough_gt           "$GROUND_TRUTH" "$FILT_GT"

# ----------------------------
# Run Truvari benches
# ----------------------------
echo "Running Truvari for ${CELL_LINE} T2T vs HG38..."
truvari bench \
  --base "$FILT_T2T" \
  --comp "$FILT_HG38" \
  --output "$HG38_T2T" \
  -f "$REFERENCE_PATH" \
  --passonly --sizemin 0

echo "Running Truvari for ${CELL_LINE} GT vs HG38..."
truvari bench \
  --base "$FILT_GT" \
  --comp "$FILT_HG38" \
  --output "$HG38_GROUNDTRUTH" \
  -f "$REFERENCE_PATH" \
  --passonly --sizemin 0

echo "Running Truvari for ${CELL_LINE} GT vs T2T..."
truvari bench \
  --base "$FILT_GT" \
  --comp "$FILT_T2T" \
  --output "$T2T_GROUNDTRUTH" \
  -f "$REFERENCE_PATH" \
  --passonly --sizemin 0

echo "Running Truvari for ${CELL_LINE} 3 groups..."
truvari bench \
  --base "$FILT_GT" \
  --comp "$HG38_T2T/tp-base.vcf.gz" \
  --output "$ALL3" \
  -f "$REFERENCE_PATH" \
  --passonly --sizemin 0

echo "Completed processing for $CELL_LINE"

conda deactivate

rm *.out
rm *.err