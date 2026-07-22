#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run PHI-base and UniProt fungi annotation using DIAMOND BLASTP
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

    [[ -f "${PHI_BASE_DMND}" ]] || die "PHI-base database not found: ${PHI_BASE_DMND}"
    [[ -f "${UNIPROT_FUNGI_DMND}" ]] || die "UniProt fungi database not found: ${UNIPROT_FUNGI_DMND}"
    [[ -f "${CM_PROTEINS}" ]] || die "CM protein file not found"
    [[ -f "${CZ_PROTEINS}" ]] || die "CZ protein file not found"

    mkdir -p "${PHI_DIR}"
}

run_phi_one() {
    local species="$1"
    local query="$2"

    log "Running PHI-base DIAMOND search for ${species}"

    diamond blastp \
        --db "${PHI_BASE_DMND}" \
        --query "${query}" \
        --out "${PHI_DIR}/${species}_phi_diamond.txt" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --evalue "${EVALUE}" \
        --max-target-seqs "${MAX_TARGET_SEQS_PHI}" \
        --threads "${THREADS_SMALL}" \
        --more-sensitive \
        > "${LOG_DIR}/${species}.phi.log" 2>&1

    log "Running UniProt fungi DIAMOND search for ${species}"

    diamond blastp \
        --db "${UNIPROT_FUNGI_DMND}" \
        --query "${query}" \
        --out "${PHI_DIR}/${species}_uniprot_fungi.txt" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --evalue "${EVALUE}" \
        --max-target-seqs "${MAX_TARGET_SEQS_UNIPROT}" \
        --threads "${THREADS_SMALL}" \
        --more-sensitive \
        > "${LOG_DIR}/${species}.uniprot_fungi.log" 2>&1
}

write_parser() {
    cat > "${PHI_DIR}/parse_phi_results.py" << 'PYEOF'
#!/usr/bin/env python3
import os
import re
import pandas as pd

COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "stitle"
]

KEYWORDS = {
    "virulence_factors": ["virulence", "pathogenicity", "toxin", "effector", "avirulence"],
    "cell_wall_degrading": ["cellulase", "pectinase", "xylanase", "cutinase", "polygalacturonase"],
    "secondary_metabolites": ["polyketide", "nonribosomal", "mycotoxin", "aflatoxin", "trichothecene"],
    "host_interaction": ["adhesin", "invasin", "penetration", "appressorium", "haustorium"],
    "stress_response": ["oxidative", "osmotic", "heat shock", "stress response", "catalase"],
    "secreted_proteins": ["secreted", "signal peptide", "extracellular", "secretory"],
    "transcription_factors": ["transcription factor", "regulator", "activator", "repressor"],
    "kinases_phosphatases": ["kinase", "phosphatase", "protein kinase", "serine threonine"]
}

def parse_blast(path):
    if not os.path.exists(path):
        return pd.DataFrame(columns=COLUMNS)
    return pd.read_csv(path, sep="\t", names=COLUMNS)

def best_hits(df):
    if df.empty:
        return df
    return df.loc[df.groupby("qseqid")["bitscore"].idxmax()].copy()

def parse_phi(prefix):
    df = parse_blast(f"{prefix}_phi_diamond.txt")
    df = df[(df["pident"] >= 40) & (df["evalue"] <= 1e-5)].copy()

    if not df.empty:
        df["phi_gene"] = df["sseqid"].astype(str).str.extract(r"PHI:(\d+)")
        df["gene_name"] = df["stitle"].astype(str).str.extract(r"^([A-Za-z0-9_.:-]+)")
        df = best_hits(df)

    df.to_csv(f"{prefix}_phi_annotation.csv", index=False)
    genes = set(df["qseqid"].dropna().unique())
    pd.DataFrame({"gene_id": sorted(genes)}).to_csv(f"{prefix}_phi_genes.txt", index=False, header=False)
    return df, genes

def parse_uniprot(prefix):
    df = parse_blast(f"{prefix}_uniprot_fungi.txt")
    df = df[(df["pident"] >= 50) & (df["evalue"] <= 1e-10)].copy()

    if not df.empty:
        df["protein_name"] = df["stitle"].astype(str).str.extract(r"([^=]+)")
        df["organism"] = df["stitle"].astype(str).str.extract(r"OS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|$)")
        df = best_hits(df)

    df.to_csv(f"{prefix}_uniprot_annotation.csv", index=False)
    genes = set(df["qseqid"].dropna().unique())
    pd.DataFrame({"gene_id": sorted(genes)}).to_csv(f"{prefix}_uniprot_genes.txt", index=False, header=False)
    return df, genes

def classify_uniprot(df, prefix):
    category_genes = {k: set() for k in KEYWORDS}

    if not df.empty:
        for _, row in df.iterrows():
            gene = row["qseqid"]
            desc = str(row["stitle"]).lower()
            for cat, kws in KEYWORDS.items():
                if any(k.lower() in desc for k in kws):
                    category_genes[cat].add(gene)

    for cat, genes in category_genes.items():
        if genes:
            pd.DataFrame({"gene_id": sorted(genes)}).to_csv(f"{prefix}_{cat}_genes.txt", index=False, header=False)

    return category_genes

def process_species(prefix):
    phi_df, phi_genes = parse_phi(prefix)
    uniprot_df, uniprot_genes = parse_uniprot(prefix)
    classify_uniprot(uniprot_df, prefix)

    high_conf = phi_genes & uniprot_genes
    all_candidates = phi_genes | uniprot_genes

    pd.DataFrame({"gene_id": sorted(high_conf)}).to_csv(f"{prefix}_high_confidence_pathogenicity.txt", index=False, header=False)
    pd.DataFrame({"gene_id": sorted(all_candidates)}).to_csv(f"{prefix}_all_pathogenicity_candidates.txt", index=False, header=False)

    print(f"{prefix} PHI genes: {len(phi_genes)}")
    print(f"{prefix} UniProt genes: {len(uniprot_genes)}")
    print(f"{prefix} high-confidence pathogenicity genes: {len(high_conf)}")
    print(f"{prefix} all pathogenicity candidates: {len(all_candidates)}")

    return high_conf, all_candidates

def main():
    cm_high, cm_all = process_species("CM")
    cz_high, cz_all = process_species("CZ")

    shared = cm_all & cz_all
    cm_specific = cm_all - cz_all
    cz_specific = cz_all - cm_all

    max_len = max(len(shared), len(cm_specific), len(cz_specific), 1)

    pd.DataFrame({
        "CM_specific": list(cm_specific) + [""] * (max_len - len(cm_specific)),
        "CZ_specific": list(cz_specific) + [""] * (max_len - len(cz_specific)),
        "Shared": list(shared) + [""] * (max_len - len(shared))
    }).to_csv("pathogenicity_comparison.csv", index=False)

    combined = cm_all | cz_all
    combined_high = cm_high | cz_high

    pd.DataFrame({"gene_id": sorted(combined)}).to_csv("combined_pathogenicity_genes.txt", index=False, header=False)
    pd.DataFrame({"gene_id": sorted(combined_high)}).to_csv("combined_high_confidence_pathogenicity.txt", index=False, header=False)

    print(f"Combined pathogenicity candidates: {len(combined)}")
    print(f"Combined high-confidence pathogenicity genes: {len(combined_high)}")

if __name__ == "__main__":
    main()
PYEOF
}

parse_results() {
    cd "${PHI_DIR}"
    python parse_phi_results.py > "${LOG_DIR}/phi.parse.log" 2>&1
}

main() {
    check_inputs
    run_phi_one "CM" "${CM_PROTEINS}"
    run_phi_one "CZ" "${CZ_PROTEINS}"
    write_parser
    parse_results
    log "PHI-base and UniProt fungi annotation completed"
}

main "$@"
