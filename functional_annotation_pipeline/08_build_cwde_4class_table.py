#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build a four-class CWDE table from dbCAN/CAZy annotation and CM/CZ orthogroup table.

Four retained plant cell wall-degrading enzyme classes:
1. Cellulase
2. Hemicellulase
3. Pectinase
4. Cutinase

If multiple CAZy_function classes are assigned to the same orthogroup entry,
a primary function is selected using the following priority:

    Cellulase > Hemicellulase > Pectinase > Cutinase

Outputs:
    dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_table.csv
    dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_table.xlsx
    dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.csv
    dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.xlsx
    dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.tsv
    dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_summary.csv

Usage:
    python 08_build_cwde_4class_table.py

Recommended usage:
    source ./config.sh
    python 08_build_cwde_4class_table.py
"""

import os
import pandas as pd
from collections import defaultdict, Counter


# ==============================================================================
# 0. Paths
# ==============================================================================

WORK_DIR = os.environ.get("WORK_DIR", os.getcwd())

OG_FILE = os.environ.get(
    "OG_CSV",
    os.path.join(WORK_DIR, "CM_CZ_OG.expanded.csv")
)

DBCAN_DIR = os.environ.get(
    "DBCAN_DIR",
    os.path.join(WORK_DIR, "dbcan_annotation")
)

os.makedirs(DBCAN_DIR, exist_ok=True)

OUT_CSV = os.path.join(DBCAN_DIR, "CM_CZ_OG_expanded_CWDE_4class_primary_table.csv")
OUT_XLSX = os.path.join(DBCAN_DIR, "CM_CZ_OG_expanded_CWDE_4class_primary_table.xlsx")

OUT_SIMPLE_CSV = os.path.join(DBCAN_DIR, "CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.csv")
OUT_SIMPLE_XLSX = os.path.join(DBCAN_DIR, "CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.xlsx")
OUT_SIMPLE_TSV = os.path.join(DBCAN_DIR, "CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.tsv")

OUT_SUMMARY = os.path.join(DBCAN_DIR, "CM_CZ_OG_expanded_CWDE_4class_primary_summary.csv")


# ==============================================================================
# 1. CWDE classification rules
# ==============================================================================

CWDE_FAMILIES = {
    "Cellulase": [
        "GH1", "GH3", "GH5", "GH6", "GH7", "GH8", "GH9",
        "GH12", "GH44", "GH45", "GH48", "GH61", "GH74", "AA9"
    ],

    "Hemicellulase": [
        "GH10", "GH11", "GH26", "GH27", "GH35", "GH36",
        "GH43", "GH51", "GH53", "GH54", "GH62", "GH67",
        "GH74", "GH95"
    ],

    "Pectinase": [
        "GH28", "GH78", "GH88", "GH105",
        "PL1", "PL2", "PL3", "PL4", "PL9", "PL10", "PL11",
        "CE8", "CE12"
    ],

    "Cutinase": [
        "CE5", "CE16"
    ]
}

FUNCTION_PRIORITY = [
    "Cellulase",
    "Hemicellulase",
    "Pectinase",
    "Cutinase"
]

SUBSTRATE_MAP = {
    # Cellulase
    "GH1": "β-glucan/Cellulose-derived oligosaccharides",
    "GH3": "β-glucan/Xylan/Arabinoxylan-derived oligosaccharides",
    "GH5": "Cellulose/Mannan/β-glucan",
    "GH6": "Cellulose",
    "GH7": "Cellulose",
    "GH8": "Cellulose/Xylan/Chitosan",
    "GH9": "Cellulose",
    "GH12": "Cellulose/Xyloglucan",
    "GH44": "Cellulose/Xyloglucan",
    "GH45": "Cellulose",
    "GH48": "Cellulose",
    "GH61": "Cellulose oxidative cleavage",
    "AA9": "Cellulose oxidative cleavage",
    "GH74": "Xyloglucan",

    # Hemicellulase
    "GH10": "Xylan",
    "GH11": "Xylan",
    "GH26": "Mannan",
    "GH27": "Galactomannan/Galactan",
    "GH35": "Galactan/Pectin side chain",
    "GH36": "Galactan/Galactomannan",
    "GH43": "Arabinoxylan/Arabinan/Xylan side chain",
    "GH51": "Arabinan/Arabinoxylan",
    "GH53": "Arabinogalactan",
    "GH54": "Arabinan",
    "GH62": "Arabinan/Arabinoxylan",
    "GH67": "Glucuronoxylan",
    "GH95": "Fucose-containing hemicellulose",

    # Pectinase
    "GH28": "Pectin/Polygalacturonan",
    "GH78": "Rhamnogalacturonan",
    "GH88": "Pectin-derived unsaturated glucuronides",
    "GH105": "Rhamnogalacturonan",
    "PL1": "Pectin/Polygalacturonan",
    "PL2": "Pectin/Polygalacturonan",
    "PL3": "Pectin/Polygalacturonan",
    "PL4": "Pectin/Rhamnogalacturonan",
    "PL9": "Pectin",
    "PL10": "Pectin",
    "PL11": "Pectin",
    "CE8": "Pectin methyl ester",
    "CE12": "Pectin acetyl ester",

    # Cutinase
    "CE5": "Cutin",
    "CE16": "Cutin/Lipid ester"
}


# ==============================================================================
# 2. Utility functions
# ==============================================================================

def clean(x):
    if pd.isna(x):
        return ""

    x = str(x).strip()

    if x in ["", "-", "NA", "nan", "NaN", "None", "none", "NULL", "null"]:
        return ""

    return x


def join_unique(values):
    vals = [clean(v) for v in values if clean(v)]
    vals = sorted(set(vals))
    return "; ".join(vals)


def infer_cazy_function(family):
    family = clean(family)

    if not family:
        return ""

    funcs = []

    for func, fams in CWDE_FAMILIES.items():
        if family in fams:
            funcs.append(func)

    return "; ".join(sorted(set(funcs)))


def infer_cazy_substrate(family):
    family = clean(family)
    return SUBSTRATE_MAP.get(family, "")


def split_function_items(x):
    x = clean(x)

    if not x:
        return []

    return [clean(i) for i in x.split(";") if clean(i)]


def is_four_class_cwde(family):
    return bool(infer_cazy_function(family))


def choose_primary_function_from_records(records):
    funcs = set()

    for r in records:
        for f in split_function_items(r.get("function", "")):
            funcs.add(f)

    for p in FUNCTION_PRIORITY:
        if p in funcs:
            return p

    return ""


def filter_records_by_primary_function(records, primary_function):
    primary_function = clean(primary_function)

    if not primary_function:
        return records

    filtered = []

    for r in records:
        funcs = split_function_items(r.get("function", ""))

        if primary_function in funcs:
            filtered.append(r)

    return filtered


# ==============================================================================
# 3. ID handling
# ==============================================================================

def normalize_for_match(gene_id):
    """
    Conservative ID normalization for matching.

    Examples:
        transcript:KAFxxx -> KAFxxx
        transcript_KAFxxx -> KAFxxx
        KAFxxx            -> KAFxxx
        KAMxxx            -> KAMxxx
    """
    gene_id = clean(gene_id)

    if not gene_id:
        return ""

    gene_id = gene_id.replace("transcript:", "")
    gene_id = gene_id.replace("transcript_", "")

    return gene_id


def display_gene_id(gene_id):
    """
    ID form used in output.

    Rules:
        transcript:xxx -> transcript:xxx
        transcript_xxx -> transcript:xxx
        KAF/PKS/PKR    -> transcript:ID
        KAM            -> unchanged
    """
    gene_id = clean(gene_id)

    if not gene_id:
        return ""

    if gene_id.startswith("transcript:"):
        return gene_id

    if gene_id.startswith("transcript_"):
        return gene_id.replace("transcript_", "transcript:", 1)

    if gene_id.startswith(("KAF", "PKS", "PKR")):
        return f"transcript:{gene_id}"

    return gene_id


def id_aliases(gene_id):
    gene_id = clean(gene_id)

    if not gene_id:
        return set()

    base = normalize_for_match(gene_id)

    aliases = {
        gene_id,
        base,
        f"transcript:{base}",
        f"transcript_{base}"
    }

    return {x for x in aliases if clean(x)}


# ==============================================================================
# 4. Load dbCAN annotation
# ==============================================================================

def add_annotation(anno_dict, gene_id, family, method, evalue="", score="", pident=""):
    family = clean(family)

    if not gene_id or not family:
        return

    if not is_four_class_cwde(family):
        return

    rec = {
        "family": family,
        "substrate": infer_cazy_substrate(family),
        "function": infer_cazy_function(family),
        "method": method,
        "evalue": clean(evalue),
        "score": clean(score),
        "pident": clean(pident)
    }

    for aid in id_aliases(gene_id):
        anno_dict[aid].append(rec)


def load_dbcan_file(path, species, method):
    anno = defaultdict(list)

    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"[WARN] File missing or empty: {path}")
        return anno

    print(f"[INFO] Loading {species} {method}: {path}")

    df = pd.read_csv(path, low_memory=False)

    if "gene_id" in df.columns:
        gene_col = "gene_id"
    elif "qseqid" in df.columns:
        gene_col = "qseqid"
    else:
        raise ValueError(
            f"{path} does not contain gene_id or qseqid column. "
            f"Columns: {df.columns.tolist()}"
        )

    if "cazy_family" not in df.columns:
        raise ValueError(
            f"{path} does not contain cazy_family column. "
            f"Columns: {df.columns.tolist()}"
        )

    total_rows = 0
    kept_rows = 0

    for _, row in df.iterrows():
        total_rows += 1

        gene_id = row.get(gene_col, "")
        family = row.get("cazy_family", "")

        if is_four_class_cwde(family):
            kept_rows += 1

        add_annotation(
            anno_dict=anno,
            gene_id=gene_id,
            family=family,
            method=method,
            evalue=row.get("evalue", ""),
            score=row.get("score", ""),
            pident=row.get("pident", "")
        )

    print(f"[INFO] {species} {method} raw records: {total_rows}")
    print(f"[INFO] {species} {method} four-class CWDE records: {kept_rows}")
    print(f"[INFO] {species} {method} ID keys: {len(anno)}")

    return anno


def merge_annotation_dicts(*dicts):
    merged = defaultdict(list)

    for d in dicts:
        for k, records in d.items():
            merged[k].extend(records)

    return merged


def load_all_dbcan():
    cm_hmm = load_dbcan_file(
        os.path.join(DBCAN_DIR, "CM_hmm_cwde.csv"),
        "CM",
        "HMM"
    )

    cm_diamond = load_dbcan_file(
        os.path.join(DBCAN_DIR, "CM_diamond_cwde.csv"),
        "CM",
        "Diamond"
    )

    cz_hmm = load_dbcan_file(
        os.path.join(DBCAN_DIR, "CZ_hmm_cwde.csv"),
        "CZ",
        "HMM"
    )

    cz_diamond = load_dbcan_file(
        os.path.join(DBCAN_DIR, "CZ_diamond_cwde.csv"),
        "CZ",
        "Diamond"
    )

    cm_anno = merge_annotation_dicts(cm_hmm, cm_diamond)
    cz_anno = merge_annotation_dicts(cz_hmm, cz_diamond)

    print(f"[INFO] CM merged ID keys: {len(cm_anno)}")
    print(f"[INFO] CZ merged ID keys: {len(cz_anno)}")

    return cm_anno, cz_anno


def get_gene_annotations(anno_dict, gene_id):
    records = []

    for aid in id_aliases(gene_id):
        records.extend(anno_dict.get(aid, []))

    seen = set()
    uniq = []

    for r in records:
        key = (
            r.get("family", ""),
            r.get("substrate", ""),
            r.get("function", ""),
            r.get("method", "")
        )

        if key not in seen:
            seen.add(key)
            uniq.append(r)

    return uniq


# ==============================================================================
# 5. Build CWDE OG table
# ==============================================================================

def build_cwde_og_table():
    if not os.path.exists(OG_FILE):
        raise FileNotFoundError(f"Orthogroup expanded file not found: {OG_FILE}")

    print(f"[INFO] Loading OG expanded table: {OG_FILE}")

    og = pd.read_csv(OG_FILE, encoding="utf-8-sig")
    og.columns = og.columns.str.strip()

    print(f"[INFO] OG expanded rows: {len(og)}")
    print(f"[INFO] OG expanded columns: {og.columns.tolist()}")

    required = ["Orthogroup", "CM_Gene.x", "CZ_Gene.x"]

    for c in required:
        if c not in og.columns:
            raise ValueError(f"OG expanded file is missing required column: {c}")

    if "copy_type" in og.columns:
        presence_col = "copy_type"
    elif "Presence_type" in og.columns:
        presence_col = "Presence_type"
    else:
        presence_col = None

    cm_anno, cz_anno = load_all_dbcan()

    records = []

    for idx, row in og.iterrows():
        if idx % 5000 == 0:
            print(f"[INFO] Integration progress: {idx}/{len(og)}")

        og_id = clean(row["Orthogroup"])
        cm_gene_raw = clean(row["CM_Gene.x"])
        cz_gene_raw = clean(row["CZ_Gene.x"])

        presence = clean(row[presence_col]) if presence_col else ""

        cm_records_raw = get_gene_annotations(cm_anno, cm_gene_raw)
        cz_records_raw = get_gene_annotations(cz_anno, cz_gene_raw)

        if not cm_records_raw and not cz_records_raw:
            continue

        all_records_raw = cm_records_raw + cz_records_raw

        primary_function = choose_primary_function_from_records(all_records_raw)

        all_records = filter_records_by_primary_function(all_records_raw, primary_function)
        cm_records = filter_records_by_primary_function(cm_records_raw, primary_function)
        cz_records = filter_records_by_primary_function(cz_records_raw, primary_function)

        families = [r["family"] for r in all_records]
        substrates = [r["substrate"] for r in all_records]
        methods = [r["method"] for r in all_records]

        rec = {
            "Orthogroup": og_id,
            "CM_Gene.x": display_gene_id(cm_gene_raw),
            "CZ_Gene.x": display_gene_id(cz_gene_raw),
            "CAZy_family": join_unique(families),
            "CAZy_substrate": join_unique(substrates),
            "CAZy_function": primary_function,
            "Presence_type": presence,

            "CM_CAZy_family": join_unique([r["family"] for r in cm_records]),
            "CZ_CAZy_family": join_unique([r["family"] for r in cz_records]),
            "CM_CAZy_function": primary_function if cm_records else "",
            "CZ_CAZy_function": primary_function if cz_records else "",
            "CM_CAZy_substrate": join_unique([r["substrate"] for r in cm_records]),
            "CZ_CAZy_substrate": join_unique([r["substrate"] for r in cz_records]),
            "dbCAN_method": join_unique(methods)
        }

        records.append(rec)

    out = pd.DataFrame(records)

    return out


# ==============================================================================
# 6. Summary
# ==============================================================================

def write_summary(out):
    stats = []

    def add(k, v):
        stats.append({"Statistic": k, "Value": v})

    add("Total_CWDE_gene_pairs", len(out))

    if len(out) == 0:
        pd.DataFrame(stats).to_csv(OUT_SUMMARY, index=False)
        print(f"[INFO] Summary file: {OUT_SUMMARY}")
        return

    add("Unique_Orthogroups", out["Orthogroup"].nunique())
    add("Unique_CM_genes", out["CM_Gene.x"].nunique())
    add("Unique_CZ_genes", out["CZ_Gene.x"].nunique())

    if "Presence_type" in out.columns:
        for k, v in out["Presence_type"].value_counts().items():
            add(f"Presence_type_{k}", v)

    fam_counter = Counter()

    for x in out["CAZy_family"]:
        for item in str(x).split(";"):
            item = item.strip()
            if item:
                fam_counter[item] += 1

    for k, v in fam_counter.most_common():
        add(f"CAZy_family_{k}", v)

    func_counter = Counter()

    for x in out["CAZy_function"]:
        item = clean(x)
        if item:
            func_counter[item] += 1

    for k, v in func_counter.most_common():
        add(f"CAZy_function_{k}", v)

    method_counter = Counter()

    for x in out["dbCAN_method"]:
        for item in str(x).split(";"):
            item = item.strip()
            if item:
                method_counter[item] += 1

    for k, v in method_counter.most_common():
        add(f"dbCAN_method_{k}", v)

    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(OUT_SUMMARY, index=False)

    print(f"[INFO] Summary file: {OUT_SUMMARY}")


# ==============================================================================
# 7. Main
# ==============================================================================

def main():
    print("============================================================")
    print("Build four-class CWDE table")
    print("Classes: Cellulase, Hemicellulase, Pectinase, Cutinase")
    print("Priority: Cellulase > Hemicellulase > Pectinase > Cutinase")
    print("============================================================")

    print(f"[INFO] WORK_DIR: {WORK_DIR}")
    print(f"[INFO] OG_FILE: {OG_FILE}")
    print(f"[INFO] DBCAN_DIR: {DBCAN_DIR}")

    out = build_cwde_og_table()

    print(f"[INFO] Output four-class CWDE rows: {len(out)}")

    simple_cols = [
        "Orthogroup",
        "CM_Gene.x",
        "CZ_Gene.x",
        "CAZy_family",
        "CAZy_substrate",
        "CAZy_function",
        "Presence_type"
    ]

    out.to_csv(OUT_CSV, index=False)
    out.to_excel(OUT_XLSX, index=False)

    out[simple_cols].to_csv(OUT_SIMPLE_CSV, index=False)
    out[simple_cols].to_excel(OUT_SIMPLE_XLSX, index=False)
    out[simple_cols].to_csv(OUT_SIMPLE_TSV, sep="\t", index=False)

    print(f"[INFO] Full CSV: {OUT_CSV}")
    print(f"[INFO] Full Excel: {OUT_XLSX}")
    print(f"[INFO] Simple CSV: {OUT_SIMPLE_CSV}")
    print(f"[INFO] Simple Excel: {OUT_SIMPLE_XLSX}")
    print(f"[INFO] Simple TSV: {OUT_SIMPLE_TSV}")

    write_summary(out)

    print("============================================================")
    print("Four-class CWDE table completed")
    print("============================================================")


if __name__ == "__main__":
    main()
