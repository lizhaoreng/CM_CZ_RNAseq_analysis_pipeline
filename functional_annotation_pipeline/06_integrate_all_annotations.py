#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integrate CM/CZ orthogroup-level functional annotations.

Integrated annotation sources:
1. Orthogroup relationship table
2. InterProScan
3. KOfamScan
4. eggNOG-mapper
5. NCBI RefSeq fungi DIAMOND BLAST
6. PHI-base / UniProt fungi pathogenicity annotation
7. MEROPS peptidase annotation
8. dbCAN / CAZy / CWDE annotation

Usage:
    python 06_integrate_all_annotations.py

Recommended usage with config.sh:
    source ./config.sh
    python 06_integrate_all_annotations.py

Main outputs:
    integrated_annotation_results/CM_CZ_OG_integrated_all_annotations.xlsx
    integrated_annotation_results/CM_CZ_OG_integrated_all_annotations.csv
    integrated_annotation_results/CWDE_orthogroups.xlsx
    integrated_annotation_results/Peptidase_orthogroups.xlsx
    integrated_annotation_results/Pathogenicity_orthogroups.xlsx
    integrated_annotation_results/integration_summary_statistics.csv
    integrated_annotation_results/integration_report.txt
"""

import os
import re
import glob
import pandas as pd
from collections import defaultdict, Counter


# ==============================================================================
# 0. Paths
# ==============================================================================

WORK_DIR = os.environ.get("WORK_DIR", os.getcwd())

OUT_DIR = os.environ.get(
    "INTEGRATED_DIR",
    os.path.join(WORK_DIR, "integrated_annotation_results")
)

OG_CSV = os.environ.get(
    "OG_CSV",
    os.path.join(WORK_DIR, "CM_CZ_OG.expanded.csv")
)

OG_XLSX = os.environ.get(
    "OG_XLSX",
    os.path.join(WORK_DIR, "Complete_CM_CZ_Ortholog_Analysis.xlsx")
)

INTERPRO_DIR = os.environ.get(
    "INTERPRO_DIR",
    os.path.join(WORK_DIR, "annotation_results", "interpro")
)

KOFAM_DIR = os.environ.get(
    "KOFAM_DIR",
    os.path.join(WORK_DIR, "annotation_results", "kofam")
)

EGGNOG_DIR = os.environ.get(
    "EGGNOG_DIR",
    os.path.join(WORK_DIR, "annotation_results", "eggnog")
)

BLAST_DIR = os.environ.get(
    "BLAST_DIR",
    os.path.join(WORK_DIR, "orthofinder_input", "blast_results")
)

PHI_DIR = os.environ.get(
    "PHI_DIR",
    os.path.join(WORK_DIR, "phi_annotation")
)

MEROPS_DIR = os.environ.get(
    "MEROPS_DIR",
    os.path.join(WORK_DIR, "merops_annotation")
)

DBCAN_DIR = os.environ.get(
    "DBCAN_DIR",
    os.path.join(WORK_DIR, "dbcan_annotation")
)

os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 1. Utility functions
# ==============================================================================

def log(msg):
    print(f"[INFO] {msg}")


def warn(msg):
    print(f"[WARN] {msg}")


def exists_nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def clean_value(x):
    if pd.isna(x):
        return ""

    x = str(x).strip()

    if x in ["", "-", "NA", "nan", "NaN", "None", "none", "NULL", "null"]:
        return ""

    return x


def clean_protein_id_for_phi(protein_id):
    """
    ID cleaning mainly used for PHI-base matching.

    Examples:
        transcript:KAFxxxx -> KAFxxxx
        transcript_KAFxxxx -> KAFxxxx
        protein:xxxx       -> xxxx
        gene:xxxx          -> xxxx
    """
    protein_id = clean_value(protein_id)

    if not protein_id:
        return ""

    for prefix in ["transcript:", "protein:", "gene:"]:
        if protein_id.startswith(prefix):
            protein_id = protein_id[len(prefix):]
            break

    if protein_id.startswith("transcript_"):
        protein_id = protein_id.replace("transcript_", "", 1)

    return protein_id


def standardize_gene_id(gene_id):
    """
    Conservative ID standardization for integration.

    Rules:
        transcript:KAFxxxx -> transcript_KAFxxxx
        transcript_KAFxxxx -> transcript_KAFxxxx
        KAFxxxx            -> transcript_KAFxxxx
        KAMxxxx            -> KAMxxxx
        Others             -> unchanged
    """
    gene_id = clean_value(gene_id)

    if not gene_id:
        return ""

    gene_id = gene_id.replace("transcript:", "transcript_")

    if gene_id.startswith("transcript_"):
        return gene_id

    if gene_id.startswith("KAF"):
        return f"transcript_{gene_id}"

    if gene_id.startswith("KAM"):
        return gene_id

    return gene_id


def id_aliases(gene_id):
    """
    Generate alternative ID forms to improve matching across annotation tools.
    """
    gene_id = clean_value(gene_id)

    if not gene_id:
        return set()

    aliases = set()

    raw = gene_id
    std = standardize_gene_id(gene_id)
    phi_clean = clean_protein_id_for_phi(gene_id)

    aliases.add(raw)
    aliases.add(std)
    aliases.add(phi_clean)

    if raw.startswith("transcript:"):
        aliases.add(raw.replace("transcript:", "transcript_"))
        aliases.add(raw.replace("transcript:", ""))

    if raw.startswith("transcript_"):
        aliases.add(raw.replace("transcript_", "", 1))
        aliases.add(raw.replace("transcript_", "transcript:", 1))

    if phi_clean.startswith("KAF"):
        aliases.add(f"transcript_{phi_clean}")
        aliases.add(f"transcript:{phi_clean}")

    return {x for x in aliases if clean_value(x)}


def add_record_with_aliases(data, gene_id, record):
    for aid in id_aliases(gene_id):
        data[aid].append(record)


def get_records_by_alias(data, gene_id):
    records = []

    for aid in id_aliases(gene_id):
        records.extend(data.get(aid, []))

    seen = set()
    unique_records = []

    for r in records:
        key = tuple(sorted((str(k), str(v)) for k, v in r.items()))
        if key not in seen:
            seen.add(key)
            unique_records.append(r)

    return unique_records


def get_one_by_alias(data, gene_id):
    for aid in id_aliases(gene_id):
        if aid in data:
            return data[aid]
    return None


def parse_gene_list(x):
    x = clean_value(x)

    if not x:
        return []

    genes = re.split(r"[,;]\s*", x)
    genes = [standardize_gene_id(g) for g in genes if clean_value(g)]

    return [g for g in genes if g]


def join_unique(values, max_items=8):
    values = [clean_value(v) for v in values if clean_value(v)]
    values = sorted(set(values))

    if not values:
        return ""

    if max_items and len(values) > max_items:
        return "; ".join(values[:max_items]) + f"; ...(+{len(values) - max_items})"

    return "; ".join(values)


def best_by_score(records, score_key="bitscore"):
    if not records:
        return None

    def get_score(r):
        try:
            return float(r.get(score_key, 0))
        except Exception:
            return 0

    return max(records, key=get_score)


# ==============================================================================
# 2. Load orthogroup table
# ==============================================================================

def classify_presence(row):
    cm = row["CM_gene_count_calc"]
    cz = row["CZ_gene_count_calc"]

    if cm == 1 and cz == 1:
        return "Single_Copy"
    if cm > 1 and cz > 1:
        return "Multi_Copy"
    if cm == 1 and cz > 1:
        return "CM_Single_CZ_Multi"
    if cm > 1 and cz == 1:
        return "CM_Multi_CZ_Single"
    if cm >= 1 and cz == 0:
        return "CM_Specific"
    if cm == 0 and cz >= 1:
        return "CZ_Specific"

    return "Unknown"


def load_og_table():
    log("Loading orthogroup table...")

    if exists_nonempty(OG_CSV):
        df = pd.read_csv(OG_CSV, encoding="utf-8-sig")
        source = OG_CSV
    elif exists_nonempty(OG_XLSX):
        df = pd.read_excel(OG_XLSX)
        source = OG_XLSX
    else:
        raise FileNotFoundError(
            "No orthogroup table found. Expected CM_CZ_OG.expanded.csv "
            "or Complete_CM_CZ_Ortholog_Analysis.xlsx"
        )

    df.columns = df.columns.str.strip()

    log(f"Orthogroup table: {source}")
    log(f"Rows: {len(df)}")
    log(f"Columns: {df.columns.tolist()}")

    cm_candidates = [
        "C_zeae_maydis", "CM_Gene", "CM_genes", "CM_Protein",
        "CM_proteins", "CM", "CM_gene", "CM_protein", "CM_Gene.x"
    ]

    cz_candidates = [
        "C_zeina", "CZ_Gene", "CZ_genes", "CZ_Protein",
        "CZ_proteins", "CZ", "CZ_gene", "CZ_protein", "CZ_Gene.x"
    ]

    cm_col = next((c for c in cm_candidates if c in df.columns), None)
    cz_col = next((c for c in cz_candidates if c in df.columns), None)

    if cm_col is None:
        for c in df.columns:
            if "CM" in c and any(k in c for k in ["Gene", "gene", "Protein", "protein"]):
                cm_col = c
                break

    if cz_col is None:
        for c in df.columns:
            if "CZ" in c and any(k in c for k in ["Gene", "gene", "Protein", "protein"]):
                cz_col = c
                break

    if cm_col is None or cz_col is None:
        raise ValueError(f"Cannot identify CM/CZ gene columns. Columns: {df.columns.tolist()}")

    if "Orthogroup" not in df.columns:
        df.insert(0, "Orthogroup", [f"OG_{i + 1:06d}" for i in range(len(df))])

    log(f"CM gene column: {cm_col}")
    log(f"CZ gene column: {cz_col}")

    df["CM_gene_list"] = df[cm_col].apply(parse_gene_list)
    df["CZ_gene_list"] = df[cz_col].apply(parse_gene_list)

    df["CM_gene_count_calc"] = df["CM_gene_list"].apply(len)
    df["CZ_gene_count_calc"] = df["CZ_gene_list"].apply(len)

    if "Presence_type" not in df.columns:
        if "copy_type" in df.columns:
            df["Presence_type"] = df["copy_type"]
        else:
            df["Presence_type"] = df.apply(classify_presence, axis=1)

    return df, cm_col, cz_col


# ==============================================================================
# 3. Load InterProScan
# ==============================================================================

def load_interpro_one(path):
    data = defaultdict(list)

    if not exists_nonempty(path):
        warn(f"InterPro file missing or empty: {path}")
        return data

    cols = [
        "protein_id", "seq_md5", "length", "analysis", "signature_id",
        "signature_desc", "start", "end", "score", "status", "date",
        "interpro_id", "interpro_desc", "go_terms", "pathways"
    ]

    df = pd.read_csv(path, sep="\t", header=None, names=cols, low_memory=False)

    for _, r in df.iterrows():
        gid = r["protein_id"]

        record = {
            "analysis": clean_value(r.get("analysis")),
            "signature_id": clean_value(r.get("signature_id")),
            "signature_desc": clean_value(r.get("signature_desc")),
            "interpro_id": clean_value(r.get("interpro_id")),
            "interpro_desc": clean_value(r.get("interpro_desc")),
            "go_terms": clean_value(r.get("go_terms"))
        }

        add_record_with_aliases(data, gid, record)

    return data


def load_interpro():
    log("Loading InterProScan annotations...")

    cm = load_interpro_one(os.path.join(INTERPRO_DIR, "CM_interpro.tsv"))
    cz = load_interpro_one(os.path.join(INTERPRO_DIR, "CZ_interpro.tsv"))

    log(f"CM InterPro keys: {len(cm)}")
    log(f"CZ InterPro keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 4. Load KOfam
# ==============================================================================

def load_kofam_one(path):
    data = defaultdict(list)

    if not exists_nonempty(path):
        warn(f"KOfam file missing or empty: {path}")
        return data

    rows = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip().replace('"', "")

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) >= 2:
                if parts[0] == "*":
                    gene = parts[1] if len(parts) > 1 else ""
                    ko = parts[2] if len(parts) > 2 else ""
                    threshold = parts[3] if len(parts) > 3 else ""
                    score = parts[4] if len(parts) > 4 else ""
                    evalue = parts[5] if len(parts) > 5 else ""
                    desc = parts[6] if len(parts) > 6 else ""
                else:
                    gene = parts[0]
                    ko = parts[1] if len(parts) > 1 else ""
                    threshold = parts[2] if len(parts) > 2 else ""
                    score = parts[3] if len(parts) > 3 else ""
                    evalue = parts[4] if len(parts) > 4 else ""
                    desc = parts[5] if len(parts) > 5 else ""

                rows.append((gene, ko, threshold, score, evalue, desc))

    for gene, ko, threshold, score, evalue, desc in rows:
        record = {
            "KO": clean_value(ko),
            "threshold": clean_value(threshold),
            "score": clean_value(score),
            "evalue": clean_value(evalue),
            "definition": clean_value(desc)
        }

        add_record_with_aliases(data, gene, record)

    return data


def load_kofam():
    log("Loading KOfam annotations...")

    cm = load_kofam_one(os.path.join(KOFAM_DIR, "CM_kofam.tsv"))
    cz = load_kofam_one(os.path.join(KOFAM_DIR, "CZ_kofam.tsv"))

    log(f"CM KOfam keys: {len(cm)}")
    log(f"CZ KOfam keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 5. Load eggNOG
# ==============================================================================

def load_eggnog_one(path):
    data = {}

    if not exists_nonempty(path):
        warn(f"eggNOG file missing or empty: {path}")
        return data

    df = pd.read_csv(path, sep="\t", comment="#", header=None, low_memory=False)

    for _, r in df.iterrows():
        raw_gid = r.iloc[0]

        record = {
            "seed_ortholog": clean_value(r.iloc[1]) if len(r) > 1 else "",
            "evalue": clean_value(r.iloc[2]) if len(r) > 2 else "",
            "score": clean_value(r.iloc[3]) if len(r) > 3 else "",
            "OGs": clean_value(r.iloc[4]) if len(r) > 4 else "",
            "COG_category": clean_value(r.iloc[6]) if len(r) > 6 else "",
            "description": clean_value(r.iloc[7]) if len(r) > 7 else "",
            "preferred_name": clean_value(r.iloc[8]) if len(r) > 8 else "",
            "GO": clean_value(r.iloc[9]) if len(r) > 9 else "",
            "EC": clean_value(r.iloc[10]) if len(r) > 10 else "",
            "KEGG_ko": clean_value(r.iloc[11]) if len(r) > 11 else "",
            "KEGG_pathway": clean_value(r.iloc[12]) if len(r) > 12 else ""
        }

        for aid in id_aliases(raw_gid):
            data[aid] = record

    return data


def load_eggnog():
    log("Loading eggNOG annotations...")

    cm = load_eggnog_one(os.path.join(EGGNOG_DIR, "CM_eggnog.emapper.annotations"))
    cz = load_eggnog_one(os.path.join(EGGNOG_DIR, "CZ_eggnog.emapper.annotations"))

    log(f"CM eggNOG keys: {len(cm)}")
    log(f"CZ eggNOG keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 6. Load NCBI fungi DIAMOND BLAST
# ==============================================================================

def load_blast_one(path):
    data = defaultdict(list)

    files = glob.glob(path)

    if not files:
        warn(f"BLAST file not found: {path}")
        return data

    path = files[0]

    if not exists_nonempty(path):
        warn(f"BLAST file empty: {path}")
        return data

    log(f"Loading BLAST file: {path}")

    df = pd.read_csv(path, sep="\t", low_memory=False)

    qcol = "Query_ID" if "Query_ID" in df.columns else df.columns[0]
    scol = "Subject_ID" if "Subject_ID" in df.columns else df.columns[1]
    pcol = "Identity%" if "Identity%" in df.columns else df.columns[2]
    eval_col = "E_value" if "E_value" in df.columns else df.columns[10]
    bit_col = "Bit_score" if "Bit_score" in df.columns else df.columns[11]
    desc_col = "Description" if "Description" in df.columns else df.columns[-1]

    for _, r in df.iterrows():
        gid = r[qcol]

        record = {
            "subject": clean_value(r[scol]),
            "identity": clean_value(r[pcol]),
            "evalue": clean_value(r[eval_col]),
            "bitscore": clean_value(r[bit_col]),
            "description": clean_value(r[desc_col])
        }

        add_record_with_aliases(data, gid, record)

    return data


def load_ncbi_blast():
    log("Loading NCBI fungi DIAMOND BLAST annotations...")

    cm = load_blast_one(os.path.join(BLAST_DIR, "CM_proteins_cleaned_blast_with_header.txt"))
    cz = load_blast_one(os.path.join(BLAST_DIR, "CZ_proteins_cleaned_blast_with_header.txt"))

    log(f"CM BLAST keys: {len(cm)}")
    log(f"CZ BLAST keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 7. Load PHI-base / UniProt pathogenicity annotations
# ==============================================================================

def load_phi_annotation_results_one(prefix):
    data = defaultdict(list)

    phi_file = os.path.join(PHI_DIR, f"{prefix}_phi_annotation.csv")
    uniprot_file = os.path.join(PHI_DIR, f"{prefix}_uniprot_annotation.csv")

    if exists_nonempty(phi_file):
        log(f"Loading {prefix} PHI-base annotation: {phi_file}")
        df = pd.read_csv(phi_file, low_memory=False)

        for _, row in df.iterrows():
            raw_id = row.get("qseqid", "")

            record = {
                "source": "PHI-base",
                "PHI_Gene_ID": clean_value(row.get("phi_gene", "")),
                "PHI_Target_ID": clean_value(row.get("sseqid", "")),
                "PHI_Identity": clean_value(row.get("pident", "")),
                "PHI_Evalue": clean_value(row.get("evalue", "")),
                "PHI_Bitscore": clean_value(row.get("bitscore", "")),
                "PHI_Description": clean_value(row.get("stitle", "")),
                "PHI_Gene_Name": clean_value(row.get("gene_name", "")),
                "UniProt_ID": "",
                "UniProt_Identity": "",
                "UniProt_Evalue": "",
                "UniProt_Bitscore": "",
                "UniProt_Protein_Name": "",
                "UniProt_Organism": "",
                "UniProt_Description": "",
                "Pathogenicity_Classifications": ""
            }

            add_record_with_aliases(data, raw_id, record)
    else:
        warn(f"{prefix} PHI-base file missing or empty: {phi_file}")

    if exists_nonempty(uniprot_file):
        log(f"Loading {prefix} UniProt annotation: {uniprot_file}")
        df = pd.read_csv(uniprot_file, low_memory=False)

        for _, row in df.iterrows():
            raw_id = row.get("qseqid", "")

            record = {
                "source": "UniProt_fungi",
                "PHI_Gene_ID": "",
                "PHI_Target_ID": "",
                "PHI_Identity": "",
                "PHI_Evalue": "",
                "PHI_Bitscore": "",
                "PHI_Description": "",
                "PHI_Gene_Name": "",
                "UniProt_ID": clean_value(row.get("sseqid", "")),
                "UniProt_Identity": clean_value(row.get("pident", "")),
                "UniProt_Evalue": clean_value(row.get("evalue", "")),
                "UniProt_Bitscore": clean_value(row.get("bitscore", "")),
                "UniProt_Protein_Name": clean_value(row.get("protein_name", "")),
                "UniProt_Organism": clean_value(row.get("organism", "")),
                "UniProt_Description": clean_value(row.get("stitle", "")),
                "Pathogenicity_Classifications": ""
            }

            add_record_with_aliases(data, raw_id, record)
    else:
        warn(f"{prefix} UniProt file missing or empty: {uniprot_file}")

    return data


def load_pathogenicity_classifications_one(prefix):
    class_data = defaultdict(set)

    classification_files = {
        "High_Confidence_Pathogenicity": [f"{prefix}_high_confidence_pathogenicity.txt"],
        "All_Pathogenicity_Candidates": [f"{prefix}_all_pathogenicity_candidates.txt"],
        "PHI_Genes": [f"{prefix}_phi_genes.txt"],
        "UniProt_Genes": [f"{prefix}_uniprot_genes.txt"],
        "Virulence_Factors": [f"{prefix}_virulence_factors_genes.txt"],
        "Cell_Wall_Degrading": [f"{prefix}_cell_wall_degrading_genes.txt"],
        "Secondary_Metabolites": [f"{prefix}_secondary_metabolites_genes.txt"],
        "Host_Interaction": [f"{prefix}_host_interaction_genes.txt"],
        "Stress_Response": [f"{prefix}_stress_response_genes.txt"],
        "Secreted_Proteins": [f"{prefix}_secreted_proteins_genes.txt"],
        "Transcription_Factors": [f"{prefix}_transcription_factors_genes.txt"],
        "Kinases_Phosphatases": [f"{prefix}_kinases_phosphatases_genes.txt"]
    }

    for category, filenames in classification_files.items():
        for filename in filenames:
            path = os.path.join(PHI_DIR, filename)

            if not exists_nonempty(path):
                continue

            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    raw_gene_id = clean_value(line)

                    if not raw_gene_id:
                        continue

                    for aid in id_aliases(raw_gene_id):
                        class_data[aid].add(category)

    return class_data


def merge_phi_annotation_and_classification(phi_data, class_data):
    merged = defaultdict(list)

    for gid, records in phi_data.items():
        for r in records:
            merged[gid].append(r)

    for gid, categories in class_data.items():
        categories_joined = join_unique(categories, max_items=None)

        record = {
            "source": "Pathogenicity_classification",
            "PHI_Gene_ID": "",
            "PHI_Target_ID": "",
            "PHI_Identity": "",
            "PHI_Evalue": "",
            "PHI_Bitscore": "",
            "PHI_Description": "",
            "PHI_Gene_Name": "",
            "UniProt_ID": "",
            "UniProt_Identity": "",
            "UniProt_Evalue": "",
            "UniProt_Bitscore": "",
            "UniProt_Protein_Name": "",
            "UniProt_Organism": "",
            "UniProt_Description": "",
            "Pathogenicity_Classifications": categories_joined
        }

        merged[gid].append(record)

    return merged


def load_phi():
    log("Loading PHI-base / UniProt pathogenicity annotations...")

    cm_phi = load_phi_annotation_results_one("CM")
    cz_phi = load_phi_annotation_results_one("CZ")

    cm_class = load_pathogenicity_classifications_one("CM")
    cz_class = load_pathogenicity_classifications_one("CZ")

    cm = merge_phi_annotation_and_classification(cm_phi, cm_class)
    cz = merge_phi_annotation_and_classification(cz_phi, cz_class)

    log(f"CM PHI/UniProt/classification keys: {len(cm)}")
    log(f"CZ PHI/UniProt/classification keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 8. Load MEROPS
# ==============================================================================

def parse_merops_description(desc):
    desc = clean_value(desc)

    mer_id = ""
    enzyme_name = ""
    organism = ""
    family = ""
    code = ""

    m = re.search(r"(MER\d+)", desc)
    if m:
        mer_id = m.group(1)

    m = re.search(r"-\s*([^(#\[]+)\s*\(([^)]+)\)", desc)
    if m:
        enzyme_name = m.group(1).strip()
        organism = m.group(2).strip()
    else:
        m = re.search(r"-\s*([^[#~]+)", desc)
        if m:
            enzyme_name = m.group(1).strip()

    m = re.search(r"\[([A-Z]\d+\.\d+[A-Z]*)\]", desc)
    if m:
        code = m.group(1)

    m = re.search(r"#([A-Z]\d+)#", desc)
    if m:
        family = m.group(1)

    return mer_id, enzyme_name, organism, family, code


def catalytic_type_from_family(family):
    family = clean_value(family)

    if not family:
        return ""

    mp = {
        "S": "serine",
        "C": "cysteine",
        "M": "metallo",
        "A": "aspartic",
        "T": "threonine",
        "G": "glutamic",
        "U": "unknown"
    }

    return mp.get(family[0].upper(), "")


def load_merops_one(path):
    data = defaultdict(list)

    if not exists_nonempty(path):
        warn(f"MEROPS file missing or empty: {path}")
        return data

    df = pd.read_csv(path, low_memory=False)

    for _, r in df.iterrows():
        gid = r.get("qseqid", "")
        desc = clean_value(r.get("stitle", ""))

        mer_id, enzyme, organism, family_from_desc, code = parse_merops_description(desc)

        family = clean_value(r.get("peptidase_family", "")) or family_from_desc
        catalytic = catalytic_type_from_family(family)

        record = {
            "family": family,
            "code": code,
            "catalytic_type": catalytic,
            "enzyme_name": enzyme,
            "organism": organism,
            "mer_id": mer_id,
            "evalue": clean_value(r.get("evalue", "")),
            "bitscore": clean_value(r.get("bitscore", "")),
            "description": desc
        }

        add_record_with_aliases(data, gid, record)

    return data


def load_merops():
    log("Loading MEROPS annotations...")

    cm = load_merops_one(os.path.join(MEROPS_DIR, "CM_merops_annotation.csv"))
    cz = load_merops_one(os.path.join(MEROPS_DIR, "CZ_merops_annotation.csv"))

    log(f"CM MEROPS keys: {len(cm)}")
    log(f"CZ MEROPS keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 9. Load dbCAN / CAZy / CWDE
# ==============================================================================

CWDE_FAMILIES = {
    "Cellulase": ["GH1", "GH3", "GH5", "GH6", "GH7", "GH8", "GH9", "GH12", "GH44", "GH45", "GH48", "GH74", "AA9"],
    "Hemicellulase": ["GH10", "GH11", "GH26", "GH27", "GH35", "GH36", "GH43", "GH51", "GH53", "GH62", "GH67", "GH95"],
    "Pectinase": ["GH28", "GH78", "GH88", "GH105", "PL1", "PL2", "PL3", "PL4", "PL9", "PL10", "PL11", "CE8"],
    "Cutinase": ["CE5", "CE16"],
    "Esterase": ["CE1", "CE2", "CE3", "CE4", "CE6", "CE8", "CE9", "CE10", "CE12", "CE15"],
    "Chitinase": ["GH18", "GH19", "GH20"]
}


def cazy_function(families):
    funcs = []

    for fam in families:
        for func, fams in CWDE_FAMILIES.items():
            if fam in fams:
                funcs.append(func)

    return join_unique(funcs, max_items=None)


def cazy_substrate(fam):
    mp = {
        "GH5": "Cellulose",
        "GH6": "Cellulose",
        "GH7": "Cellulose",
        "GH12": "Cellulose",
        "GH45": "Cellulose",
        "AA9": "Cellulose oxidative cleavage",
        "GH10": "Xylan",
        "GH11": "Xylan",
        "GH43": "Arabinoxylan",
        "GH51": "Arabinan",
        "GH28": "Pectin",
        "PL1": "Pectin",
        "PL3": "Pectin",
        "CE8": "Pectin ester",
        "CE5": "Cutin",
        "CE16": "Cutin/Lipid ester",
        "GH18": "Chitin",
        "GH19": "Chitin",
        "GH20": "Chitin"
    }

    return mp.get(fam, "")


def load_dbcan_one(prefix):
    data = defaultdict(list)

    files = [
        os.path.join(DBCAN_DIR, f"{prefix}_hmm_cwde.csv"),
        os.path.join(DBCAN_DIR, f"{prefix}_diamond_cwde.csv")
    ]

    for path in files:
        if not exists_nonempty(path):
            warn(f"dbCAN file missing or empty: {path}")
            continue

        method = "HMM" if "hmm" in os.path.basename(path).lower() else "Diamond"
        df = pd.read_csv(path, low_memory=False)

        gene_col = "gene_id" if "gene_id" in df.columns else "qseqid"

        for _, r in df.iterrows():
            gid = r.get(gene_col, "")
            fam = clean_value(r.get("cazy_family", ""))

            if not clean_value(gid) or not fam:
                continue

            record = {
                "method": method,
                "family": fam,
                "function": cazy_function([fam]),
                "substrate": cazy_substrate(fam),
                "evalue": clean_value(r.get("evalue", "")),
                "score": clean_value(r.get("score", "")),
                "pident": clean_value(r.get("pident", ""))
            }

            add_record_with_aliases(data, gid, record)

    return data


def load_dbcan():
    log("Loading dbCAN / CAZy / CWDE annotations...")

    cm = load_dbcan_one("CM")
    cz = load_dbcan_one("CZ")

    log(f"CM dbCAN keys: {len(cm)}")
    log(f"CZ dbCAN keys: {len(cz)}")

    return cm, cz


# ==============================================================================
# 10. Summarize gene annotations
# ==============================================================================

def summarize_gene_annotations(genes, datasets):
    interpro, kofam, eggnog, blast, phi, merops, dbcan = datasets

    out = {}

    # InterPro
    ipr_ids, ipr_descs, ipr_go = [], [], []

    for g in genes:
        for r in get_records_by_alias(interpro, g):
            ipr_ids.append(r.get("interpro_id", ""))
            ipr_descs.append(r.get("interpro_desc", "") or r.get("signature_desc", ""))
            ipr_go.append(r.get("go_terms", ""))

    out["InterPro_IDs"] = join_unique(ipr_ids)
    out["InterPro_desc"] = join_unique(ipr_descs)
    out["InterPro_GO"] = join_unique(ipr_go)

    # KOfam
    kos, ko_descs = [], []

    for g in genes:
        records = get_records_by_alias(kofam, g)
        best = best_by_score(records, "score")

        if best:
            kos.append(best.get("KO", ""))
            ko_descs.append(best.get("definition", ""))

    out["KO"] = join_unique(kos)
    out["KO_desc"] = join_unique(ko_descs)

    # eggNOG
    egg_desc, egg_go, egg_ec, egg_ko, egg_pathway, cog = [], [], [], [], [], []

    for g in genes:
        r = get_one_by_alias(eggnog, g)

        if r:
            egg_desc.append(r.get("description", ""))
            egg_go.append(r.get("GO", ""))
            egg_ec.append(r.get("EC", ""))
            egg_ko.append(r.get("KEGG_ko", ""))
            egg_pathway.append(r.get("KEGG_pathway", ""))
            cog.append(r.get("COG_category", ""))

    out["eggNOG_desc"] = join_unique(egg_desc)
    out["eggNOG_GO"] = join_unique(egg_go)
    out["eggNOG_EC"] = join_unique(egg_ec)
    out["eggNOG_KEGG_ko"] = join_unique(egg_ko)
    out["eggNOG_pathway"] = join_unique(egg_pathway)
    out["COG_category"] = join_unique(cog)

    # NCBI BLAST
    blast_desc, blast_subject, blast_identity, blast_evalue = [], [], [], []

    for g in genes:
        records = get_records_by_alias(blast, g)
        best = best_by_score(records, "bitscore")

        if best:
            blast_subject.append(best.get("subject", ""))
            blast_desc.append(best.get("description", ""))
            blast_identity.append(best.get("identity", ""))
            blast_evalue.append(best.get("evalue", ""))

    out["NCBI_best_hit"] = join_unique(blast_subject)
    out["NCBI_desc"] = join_unique(blast_desc)
    out["NCBI_identity"] = join_unique(blast_identity)
    out["NCBI_evalue"] = join_unique(blast_evalue)

    # PHI / UniProt
    phi_sources, phi_gene_ids, phi_target_ids, phi_gene_names = [], [], [], []
    phi_descs, phi_identity, phi_evalue = [], [], []

    uniprot_ids, uniprot_names, uniprot_orgs = [], [], []
    uniprot_descs, uniprot_identity, uniprot_evalue = [], [], []

    patho_classes = []

    for g in genes:
        records = get_records_by_alias(phi, g)

        for r in records:
            phi_sources.append(r.get("source", ""))

            phi_gene_ids.append(r.get("PHI_Gene_ID", ""))
            phi_target_ids.append(r.get("PHI_Target_ID", ""))
            phi_gene_names.append(r.get("PHI_Gene_Name", ""))
            phi_descs.append(r.get("PHI_Description", ""))
            phi_identity.append(r.get("PHI_Identity", ""))
            phi_evalue.append(r.get("PHI_Evalue", ""))

            uniprot_ids.append(r.get("UniProt_ID", ""))
            uniprot_names.append(r.get("UniProt_Protein_Name", ""))
            uniprot_orgs.append(r.get("UniProt_Organism", ""))
            uniprot_descs.append(r.get("UniProt_Description", ""))
            uniprot_identity.append(r.get("UniProt_Identity", ""))
            uniprot_evalue.append(r.get("UniProt_Evalue", ""))

            if r.get("Pathogenicity_Classifications"):
                for item in str(r.get("Pathogenicity_Classifications", "")).split(";"):
                    item = item.strip()
                    if item:
                        patho_classes.append(item)

    out["Pathogenicity_source"] = join_unique(phi_sources, max_items=None)
    out["PHI_gene"] = join_unique(phi_gene_ids, max_items=None)
    out["PHI_target"] = join_unique(phi_target_ids)
    out["PHI_gene_name"] = join_unique(phi_gene_names)
    out["PHI_identity"] = join_unique(phi_identity)
    out["PHI_evalue"] = join_unique(phi_evalue)
    out["PHI_desc"] = join_unique(phi_descs)

    out["UniProt_ID"] = join_unique(uniprot_ids)
    out["UniProt_protein_name"] = join_unique(uniprot_names)
    out["UniProt_organism"] = join_unique(uniprot_orgs)
    out["UniProt_identity"] = join_unique(uniprot_identity)
    out["UniProt_evalue"] = join_unique(uniprot_evalue)
    out["UniProt_desc"] = join_unique(uniprot_descs)

    out["Pathogenicity_classification"] = join_unique(patho_classes, max_items=None)
    out["Pathogenicity_desc"] = join_unique(phi_descs + uniprot_descs + patho_classes, max_items=12)

    out["Has_pathogenicity_hit"] = bool(
        phi_gene_ids or phi_target_ids or phi_descs or
        uniprot_ids or uniprot_descs or patho_classes
    )

    # MEROPS
    mer_fam, mer_cat, mer_enzyme, mer_id = [], [], [], []

    for g in genes:
        records = get_records_by_alias(merops, g)
        best = best_by_score(records, "bitscore")

        if best:
            mer_fam.append(best.get("family", ""))
            mer_cat.append(best.get("catalytic_type", ""))
            mer_enzyme.append(best.get("enzyme_name", ""))
            mer_id.append(best.get("mer_id", ""))

    out["MEROPS_family"] = join_unique(mer_fam)
    out["MEROPS_catalytic_type"] = join_unique(mer_cat)
    out["MEROPS_enzyme"] = join_unique(mer_enzyme)
    out["MEROPS_ID"] = join_unique(mer_id)
    out["Has_peptidase"] = bool(mer_fam or mer_cat or mer_enzyme)

    # dbCAN / CAZy / CWDE
    cazy_fam, cazy_func, cazy_sub, dbcan_methods = [], [], [], []

    for g in genes:
        records = get_records_by_alias(dbcan, g)

        for r in records:
            cazy_fam.append(r.get("family", ""))
            cazy_func.append(r.get("function", ""))
            cazy_sub.append(r.get("substrate", ""))
            dbcan_methods.append(r.get("method", ""))

    out["CAZy_family"] = join_unique(cazy_fam, max_items=None)
    out["CAZy_function"] = join_unique(cazy_func, max_items=None)
    out["CAZy_substrate"] = join_unique(cazy_sub, max_items=None)
    out["dbCAN_methods"] = join_unique(dbcan_methods, max_items=None)
    out["Has_CWDE"] = bool(cazy_fam)

    return out


# ==============================================================================
# 11. Integrate annotations
# ==============================================================================

def integrate():
    og_df, cm_col, cz_col = load_og_table()

    cm_interpro, cz_interpro = load_interpro()
    cm_kofam, cz_kofam = load_kofam()
    cm_eggnog, cz_eggnog = load_eggnog()
    cm_blast, cz_blast = load_ncbi_blast()
    cm_phi, cz_phi = load_phi()
    cm_merops, cz_merops = load_merops()
    cm_dbcan, cz_dbcan = load_dbcan()

    cm_datasets = (
        cm_interpro, cm_kofam, cm_eggnog, cm_blast,
        cm_phi, cm_merops, cm_dbcan
    )

    cz_datasets = (
        cz_interpro, cz_kofam, cz_eggnog, cz_blast,
        cz_phi, cz_merops, cz_dbcan
    )

    log("Integrating annotations at orthogroup level...")

    records = []

    for i, row in og_df.iterrows():
        if i % 1000 == 0:
            log(f"Progress: {i}/{len(og_df)}")

        cm_genes = row["CM_gene_list"]
        cz_genes = row["CZ_gene_list"]

        cm_sum = summarize_gene_annotations(cm_genes, cm_datasets)
        cz_sum = summarize_gene_annotations(cz_genes, cz_datasets)

        rec = row.drop(labels=["CM_gene_list", "CZ_gene_list"]).to_dict()

        rec["CM_genes_standardized"] = join_unique(cm_genes, max_items=None)
        rec["CZ_genes_standardized"] = join_unique(cz_genes, max_items=None)

        for k, v in cm_sum.items():
            rec[f"CM_{k}"] = v

        for k, v in cz_sum.items():
            rec[f"CZ_{k}"] = v

        combined_keys = [
            "InterPro_IDs", "InterPro_desc", "InterPro_GO",
            "KO", "KO_desc",
            "eggNOG_desc", "eggNOG_GO", "eggNOG_EC",
            "eggNOG_KEGG_ko", "eggNOG_pathway", "COG_category",
            "NCBI_best_hit", "NCBI_desc", "NCBI_identity", "NCBI_evalue",
            "Pathogenicity_source",
            "PHI_gene", "PHI_target", "PHI_gene_name", "PHI_identity",
            "PHI_evalue", "PHI_desc",
            "UniProt_ID", "UniProt_protein_name", "UniProt_organism",
            "UniProt_identity", "UniProt_evalue", "UniProt_desc",
            "Pathogenicity_classification", "Pathogenicity_desc",
            "MEROPS_family", "MEROPS_catalytic_type", "MEROPS_enzyme", "MEROPS_ID",
            "CAZy_family", "CAZy_function", "CAZy_substrate", "dbCAN_methods"
        ]

        for k in combined_keys:
            rec[f"Combined_{k}"] = join_unique(
                [cm_sum.get(k, ""), cz_sum.get(k, "")],
                max_items=None
            )

        rec["Combined_Has_pathogenicity_hit"] = bool(
            cm_sum["Has_pathogenicity_hit"] or cz_sum["Has_pathogenicity_hit"]
        )

        rec["Combined_Has_peptidase"] = bool(
            cm_sum["Has_peptidase"] or cz_sum["Has_peptidase"]
        )

        rec["Combined_Has_CWDE"] = bool(
            cm_sum["Has_CWDE"] or cz_sum["Has_CWDE"]
        )

        tags = []

        if rec["Combined_Has_CWDE"]:
            tags.append("CWDE")
        if rec["Combined_Has_peptidase"]:
            tags.append("Peptidase")
        if rec["Combined_Has_pathogenicity_hit"]:
            tags.append("Pathogenicity")
        if rec["Combined_KO"] or rec["Combined_eggNOG_KEGG_ko"]:
            tags.append("KEGG")
        if rec["Combined_InterPro_IDs"]:
            tags.append("InterPro")
        if rec["Combined_NCBI_desc"]:
            tags.append("NCBI_BLAST")

        rec["Functional_tags"] = "; ".join(tags)

        score = 0

        if rec["Combined_InterPro_IDs"]:
            score += 2
        if rec["Combined_KO"]:
            score += 2
        if rec["Combined_eggNOG_desc"]:
            score += 2
        if rec["Combined_NCBI_desc"]:
            score += 1
        if rec["Combined_Has_CWDE"]:
            score += 3
        if rec["Combined_Has_peptidase"]:
            score += 3
        if rec["Combined_Has_pathogenicity_hit"]:
            score += 3

        rec["Annotation_quality_score"] = score

        records.append(rec)

    integrated = pd.DataFrame(records)

    return integrated


# ==============================================================================
# 12. Write outputs
# ==============================================================================

def write_outputs(df):
    log("Writing integrated annotation outputs...")

    main_xlsx = os.path.join(OUT_DIR, "CM_CZ_OG_integrated_all_annotations.xlsx")
    main_csv = os.path.join(OUT_DIR, "CM_CZ_OG_integrated_all_annotations.csv")

    df.to_excel(main_xlsx, index=False)
    df.to_csv(main_csv, index=False)

    log(f"Main Excel: {main_xlsx}")
    log(f"Main CSV: {main_csv}")

    subsets = {
        "CWDE": df[df["Combined_Has_CWDE"] == True].copy(),
        "Peptidase": df[df["Combined_Has_peptidase"] == True].copy(),
        "Pathogenicity": df[df["Combined_Has_pathogenicity_hit"] == True].copy(),
        "High_quality": df[df["Annotation_quality_score"] >= 8].copy(),
        "CM_specific": df[df["Presence_type"] == "CM_Specific"].copy() if "Presence_type" in df.columns else pd.DataFrame(),
        "CZ_specific": df[df["Presence_type"] == "CZ_Specific"].copy() if "Presence_type" in df.columns else pd.DataFrame(),
    }

    for name, sub in subsets.items():
        if len(sub) > 0:
            path = os.path.join(OUT_DIR, f"{name}_orthogroups.xlsx")
            sub.to_excel(path, index=False)
            log(f"{name} table: {path} ({len(sub)} rows)")

    if "Combined_Pathogenicity_classification" in df.columns:
        high_phi = df[
            df["Combined_Pathogenicity_classification"]
            .fillna("")
            .astype(str)
            .str.contains("High_Confidence_Pathogenicity", na=False)
        ].copy()

        if len(high_phi) > 0:
            path = os.path.join(OUT_DIR, "High_confidence_pathogenicity_orthogroups.xlsx")
            high_phi.to_excel(path, index=False)
            log(f"High-confidence pathogenicity table: {path} ({len(high_phi)} rows)")

    stats = []

    def add_stat(k, v):
        stats.append({"Statistic": k, "Value": v})

    add_stat("Total_Orthogroups", len(df))

    if "Presence_type" in df.columns:
        for k, v in df["Presence_type"].value_counts().items():
            add_stat(f"Presence_type_{k}", v)

    add_stat("With_InterPro", (df["Combined_InterPro_IDs"] != "").sum())
    add_stat("With_KOfam", (df["Combined_KO"] != "").sum())
    add_stat("With_eggNOG", (df["Combined_eggNOG_desc"] != "").sum())
    add_stat("With_NCBI_BLAST", (df["Combined_NCBI_desc"] != "").sum())

    add_stat("With_PHI_base", (df["Combined_PHI_desc"] != "").sum())
    add_stat("With_UniProt_pathogenicity", (df["Combined_UniProt_desc"] != "").sum())
    add_stat("With_Pathogenicity_classification", (df["Combined_Pathogenicity_classification"] != "").sum())

    add_stat("With_CWDE", df["Combined_Has_CWDE"].sum())
    add_stat("With_Peptidase", df["Combined_Has_peptidase"].sum())
    add_stat("With_Pathogenicity", df["Combined_Has_pathogenicity_hit"].sum())
    add_stat("High_quality_score_ge_8", (df["Annotation_quality_score"] >= 8).sum())

    if "Combined_Pathogenicity_classification" in df.columns:
        patho_counter = Counter()

        for x in df["Combined_Pathogenicity_classification"]:
            for item in str(x).split(";"):
                item = item.strip()
                if item:
                    patho_counter[item] += 1

        for k, v in patho_counter.most_common():
            add_stat(f"Pathogenicity_class_{k}", v)

    cazy_counter = Counter()

    for x in df["Combined_CAZy_function"]:
        for item in str(x).split(";"):
            item = item.strip()
            if item:
                cazy_counter[item] += 1

    for k, v in cazy_counter.most_common():
        add_stat(f"CAZy_function_{k}", v)

    mer_counter = Counter()

    for x in df["Combined_MEROPS_catalytic_type"]:
        for item in str(x).split(";"):
            item = item.strip()
            if item:
                mer_counter[item] += 1

    for k, v in mer_counter.most_common():
        add_stat(f"MEROPS_catalytic_{k}", v)

    stats_df = pd.DataFrame(stats)
    stats_path = os.path.join(OUT_DIR, "integration_summary_statistics.csv")
    stats_df.to_csv(stats_path, index=False)

    log(f"Summary statistics: {stats_path}")

    report_path = os.path.join(OUT_DIR, "integration_report.txt")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("CM/CZ orthogroup functional annotation integration report\n")
        f.write("=========================================================\n\n")

        f.write(f"WORK_DIR: {WORK_DIR}\n")
        f.write(f"OG_CSV: {OG_CSV}\n")
        f.write(f"OG_XLSX: {OG_XLSX}\n")
        f.write(f"INTERPRO_DIR: {INTERPRO_DIR}\n")
        f.write(f"KOFAM_DIR: {KOFAM_DIR}\n")
        f.write(f"EGGNOG_DIR: {EGGNOG_DIR}\n")
        f.write(f"BLAST_DIR: {BLAST_DIR}\n")
        f.write(f"PHI_DIR: {PHI_DIR}\n")
        f.write(f"MEROPS_DIR: {MEROPS_DIR}\n")
        f.write(f"DBCAN_DIR: {DBCAN_DIR}\n")
        f.write(f"OUT_DIR: {OUT_DIR}\n\n")

        for item in stats:
            f.write(f"{item['Statistic']}: {item['Value']}\n")

    log(f"Report: {report_path}")


# ==============================================================================
# 13. Main
# ==============================================================================

def main():
    log("Starting full annotation integration...")
    df = integrate()
    write_outputs(df)
    log("All annotation integration completed.")


if __name__ == "__main__":
    main()
