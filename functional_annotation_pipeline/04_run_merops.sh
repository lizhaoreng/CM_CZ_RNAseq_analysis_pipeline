#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run MEROPS peptidase annotation using DIAMOND BLASTP
# ==============================================================================

source ./config.sh

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_inputs() {
    command -v diamond >/dev/null 2>&1 || die "diamond not found"
    command -v python >/dev/null 2>&1 || die "python not found"

    [[ -f "${MEROPS_DMND}" ]] || die "MEROPS database not found: ${MEROPS_DMND}"
    [[ -f "${CM_PROTEINS}" ]] || die "CM protein file not found"
    [[ -f "${CZ_PROTEINS}" ]] || die "CZ protein file not found"

    mkdir -p "${MEROPS_DIR}"
}

run_merops_one() {
    local species="$1"
    local query="$2"
    local output="${MEROPS_DIR}/${species}_merops_diamond.txt"
    local log_file="${LOG_DIR}/${species}.merops.log"

    log "Running MEROPS DIAMOND BLAST for ${species}"

    diamond blastp \
        --db "${MEROPS_DMND}" \
        --query "${query}" \
        --out "${output}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --evalue "${EVALUE}" \
        --max-target-seqs "${MAX_TARGET_SEQS_MEROPS}" \
        --threads "${THREADS_SMALL}" \
        --more-sensitive \
        > "${log_file}" 2>&1

    log "${species} MEROPS hits: $(wc -l < "${output}")"
}

write_parser() {
    cat > "${MEROPS_DIR}/parse_merops_results.py" << 'PYEOF'
#!/usr/bin/env python3
import os
import re
import pandas as pd

COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "stitle"
]

def parse_merops_blast(blast_file, species):
    if not os.path.exists(blast_file):
        print(f"Missing file: {blast_file}")
        return pd.DataFrame()

    df = pd.read_csv(blast_file, sep="\t", names=COLUMNS)
    print(f"{species}: raw hits = {len(df)}")

    df = df[(df["pident"] >= 30) & (df["evalue"] <= 1e-5)].copy()
    print(f"{species}: filtered hits = {len(df)}")

    if df.empty:
        return df

    df["peptidase_family"] = df["stitle"].astype(str).str.extract(r"(M\d+|S\d+|C\d+|A\d+|T\d+|G\d+|P\d+|U\d+)")
    df["peptidase_type"] = df["stitle"].astype(str).str.extract(
        r"(metallopeptidase|serine peptidase|cysteine peptidase|aspartic peptidase|threonine peptidase|glutamic peptidase)",
        flags=re.IGNORECASE
    )

    best_hits = df.loc[df.groupby("qseqid")["bitscore"].idxmax()].copy()
    return best_hits

def extract_peptidases(df, prefix):
    if df.empty:
        return set()

    genes = set(df["qseqid"].dropna().unique())

    df.to_csv(f"{prefix}_merops_annotation.csv", index=False)
    pd.DataFrame({"gene_id": sorted(genes)}).to_csv(f"{prefix}_peptidase_genes.txt", index=False, header=False)

    family_stats = df["peptidase_family"].value_counts(dropna=True)
    family_stats.to_csv(f"{prefix}_merops_family_summary.csv")

    return genes

def main():
    cm_df = parse_merops_blast("CM_merops_diamond.txt", "CM")
    cz_df = parse_merops_blast("CZ_merops_diamond.txt", "CZ")

    cm_genes = extract_peptidases(cm_df, "CM")
    cz_genes = extract_peptidases(cz_df, "CZ")

    shared = cm_genes & cz_genes
    cm_specific = cm_genes - cz_genes
    cz_specific = cz_genes - cm_genes

    max_len = max(len(shared), len(cm_specific), len(cz_specific), 1)

    pd.DataFrame({
        "CM_specific": list(cm_specific) + [""] * (max_len - len(cm_specific)),
        "CZ_specific": list(cz_specific) + [""] * (max_len - len(cz_specific)),
        "Shared": list(shared) + [""] * (max_len - len(shared))
    }).to_csv("peptidase_comparison.csv", index=False)

    all_genes = cm_genes | cz_genes
    pd.DataFrame({"gene_id": sorted(all_genes)}).to_csv("combined_peptidase_genes.txt", index=False, header=False)

    print(f"CM peptidases: {len(cm_genes)}")
    print(f"CZ peptidases: {len(cz_genes)}")
    print(f"Combined peptidases: {len(all_genes)}")

if __name__ == "__main__":
    main()
PYEOF
}

parse_results() {
    cd "${MEROPS_DIR}"
    python parse_merops_results.py > "${LOG_DIR}/merops.parse.log" 2>&1
}

main() {
    check_inputs
    run_merops_one "CM" "${CM_PROTEINS}"
    run_merops_one "CZ" "${CZ_PROTEINS}"
    write_parser
    parse_results
    log "MEROPS annotation completed"
}

main "$@"
