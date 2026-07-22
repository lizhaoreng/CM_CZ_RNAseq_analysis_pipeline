#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run dbCAN-like CAZy annotation using HMMER and DIAMOND, then parse CWDE genes
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
    command -v hmmscan >/dev/null 2>&1 || die "hmmscan not found"
    command -v hmmpress >/dev/null 2>&1 || die "hmmpress not found"
    command -v diamond >/dev/null 2>&1 || die "diamond not found"
    command -v python >/dev/null 2>&1 || die "python not found"

    [[ -f "${CM_PROTEINS}" ]] || die "CM protein file not found: ${CM_PROTEINS}"
    [[ -f "${CZ_PROTEINS}" ]] || die "CZ protein file not found: ${CZ_PROTEINS}"
    [[ -f "${DBCAN_HMM_DB}" ]] || die "dbCAN HMM database not found: ${DBCAN_HMM_DB}"
    [[ -f "${CAZY_DB_FA}" ]] || die "CAZy database FASTA not found: ${CAZY_DB_FA}"

    mkdir -p "${DBCAN_DIR}"
}

prepare_dbcan_db() {
    log "Preparing dbCAN databases..."

    if [[ ! -f "${DBCAN_HMM_DB}.h3m" ]]; then
        log "Running hmmpress for dbCAN HMM database..."
        hmmpress "${DBCAN_HMM_DB}"
    fi

    if [[ ! -f "${CAZY_DB_DMND}" ]]; then
        log "Building DIAMOND database for CAZy..."
        diamond makedb --in "${CAZY_DB_FA}" -d "${CAZY_DB_DMND}"
    fi
}

run_dbcan_one() {
    local species="$1"
    local query="$2"
    local outdir="${DBCAN_DIR}/${species}_dbcan_output"

    mkdir -p "${outdir}"

    log "Running dbCAN HMM search for ${species}"
    hmmscan \
        --domtblout "${outdir}/hmmscan.out" \
        --cpu "${THREADS_SMALL}" \
        "${DBCAN_HMM_DB}" \
        "${query}" \
        > "${LOG_DIR}/${species}.dbcan.hmmscan.log" 2>&1

    log "Running CAZy DIAMOND search for ${species}"
    diamond blastp \
        --db "${CAZY_DB_DMND}" \
        --query "${query}" \
        --out "${outdir}/diamond.out" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue "${EVALUE}" \
        --max-target-seqs 10 \
        --threads "${THREADS_SMALL}" \
        > "${LOG_DIR}/${species}.dbcan.diamond.log" 2>&1

    log "${species} dbCAN/CAZy search completed"
}

write_parser() {
    cat > "${DBCAN_DIR}/parse_dbcan_results.py" << 'PYEOF'
#!/usr/bin/env python3
import os
import re
import pandas as pd

CWDE_FAMILIES = {
    "Cellulase": ["GH1", "GH3", "GH5", "GH6", "GH7", "GH8", "GH9", "GH12", "GH44", "GH45", "GH48", "GH61", "GH74", "AA9"],
    "Hemicellulase": ["GH10", "GH11", "GH26", "GH27", "GH35", "GH36", "GH43", "GH51", "GH53", "GH54", "GH62", "GH67", "GH95"],
    "Pectinase": ["GH28", "GH78", "GH88", "GH105", "PL1", "PL2", "PL3", "PL4", "PL9", "PL10", "PL11", "CE8", "CE12"],
    "Cutinase": ["CE5", "CE16"],
    "Esterase": ["CE1", "CE2", "CE3", "CE4", "CE6", "CE8", "CE9", "CE10", "CE12", "CE15"],
    "Chitinase": ["GH18", "GH19", "GH20"]
}

ALL_CWDE_FAMILIES = sorted({f for fams in CWDE_FAMILIES.values() for f in fams})

def parse_hmmscan_output(hmm_file):
    results = []
    if not os.path.exists(hmm_file):
        return pd.DataFrame()

    with open(hmm_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 23:
                target_name = parts[0]
                query_name = parts[3]
                evalue = float(parts[6])
                score = float(parts[7])
                bias = float(parts[8])
                cazy_family = target_name.split(".")[0]
                results.append({
                    "gene_id": query_name,
                    "cazy_family": cazy_family,
                    "target_name": target_name,
                    "evalue": evalue,
                    "score": score,
                    "bias": bias
                })
    return pd.DataFrame(results)

def parse_diamond_output(diamond_file):
    if not os.path.exists(diamond_file):
        return pd.DataFrame()

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    df = pd.read_csv(diamond_file, sep="\t", names=columns)
    df["cazy_family"] = df["sseqid"].astype(str).str.extract(r"(GH\d+|GT\d+|PL\d+|CE\d+|AA\d+|CBM\d+)")
    return df

def extract_cwde_genes(hmm_df, diamond_df, prefix):
    hmm_cwde = hmm_df[hmm_df["cazy_family"].isin(ALL_CWDE_FAMILIES)].copy()
    diamond_cwde = diamond_df[diamond_df["cazy_family"].isin(ALL_CWDE_FAMILIES)].copy()

    hmm_filtered = hmm_cwde[(hmm_cwde["evalue"] <= 1e-5) & (hmm_cwde["score"] >= 30)].copy()
    diamond_filtered = diamond_cwde[(diamond_cwde["evalue"] <= 1e-5) & (diamond_cwde["pident"] >= 30)].copy()

    hmm_genes = set(hmm_filtered["gene_id"].dropna().unique())
    diamond_genes = set(diamond_filtered["qseqid"].dropna().unique())

    high_confidence = hmm_genes & diamond_genes
    all_candidates = hmm_genes | diamond_genes

    print(f"{prefix} HMM CWDE genes: {len(hmm_genes)}")
    print(f"{prefix} DIAMOND CWDE genes: {len(diamond_genes)}")
    print(f"{prefix} high-confidence CWDE genes: {len(high_confidence)}")
    print(f"{prefix} all candidate CWDE genes: {len(all_candidates)}")

    hmm_filtered.to_csv(f"{prefix}_hmm_cwde.csv", index=False)
    diamond_filtered.to_csv(f"{prefix}_diamond_cwde.csv", index=False)

    pd.DataFrame({"gene_id": sorted(high_confidence)}).to_csv(f"{prefix}_cwde_high_confidence.txt", index=False, header=False)
    pd.DataFrame({"gene_id": sorted(all_candidates)}).to_csv(f"{prefix}_cwde_all_candidates.txt", index=False, header=False)

    return high_confidence, all_candidates

def main():
    cm_hmm = parse_hmmscan_output("CM_dbcan_output/hmmscan.out")
    cm_diamond = parse_diamond_output("CM_dbcan_output/diamond.out")
    cz_hmm = parse_hmmscan_output("CZ_dbcan_output/hmmscan.out")
    cz_diamond = parse_diamond_output("CZ_dbcan_output/diamond.out")

    cm_high, cm_all = extract_cwde_genes(cm_hmm, cm_diamond, "CM")
    cz_high, cz_all = extract_cwde_genes(cz_hmm, cz_diamond, "CZ")

    shared = cm_high & cz_high
    cm_specific = cm_high - cz_high
    cz_specific = cz_high - cm_high

    max_len = max(len(shared), len(cm_specific), len(cz_specific), 1)

    comparison = pd.DataFrame({
        "CM_specific": list(cm_specific) + [""] * (max_len - len(cm_specific)),
        "CZ_specific": list(cz_specific) + [""] * (max_len - len(cz_specific)),
        "Shared": list(shared) + [""] * (max_len - len(shared))
    })

    comparison.to_csv("cwde_comparison.csv", index=False)

    all_cwde = cm_high | cz_high
    pd.DataFrame({"gene_id": sorted(all_cwde)}).to_csv("combined_cwde_genes.txt", index=False, header=False)

    print(f"Combined high-confidence CWDE genes: {len(all_cwde)}")

if __name__ == "__main__":
    main()
PYEOF
}

parse_results() {
    log "Parsing dbCAN results..."
    cd "${DBCAN_DIR}"
    python parse_dbcan_results.py > "${LOG_DIR}/dbcan.parse.log" 2>&1
}

main() {
    check_inputs
    prepare_dbcan_db

    run_dbcan_one "CM" "${CM_PROTEINS}"
    run_dbcan_one "CZ" "${CZ_PROTEINS}"

    write_parser
    parse_results

    log "dbCAN/CWDE annotation completed"
    log "Output directory: ${DBCAN_DIR}"
}

main "$@"
