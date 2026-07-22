#!/usr/bin/env python3

"""
01_prepare_ortholog_data.py

Prepare one-to-one orthologous gene pairs for dN/dS analysis.

This script:
  1. Reads a one-to-one ortholog table.
  2. Normalizes protein/CDS identifiers from different FASTA header formats.
  3. Generates a mapping table from original IDs to standardized orthogroup IDs.
  4. Renames CDS and protein FASTA records for downstream pair extraction.

Supported FASTA ID formats include:
  - KAF2216363
  - transcript:KAF2216363
  - KAM3422873.1
  - lcl|MVDW03000002.1_cds_KAM3422499.1_1
  - cds-KAM3422499.1
"""

import argparse
import os
import re
import sys
import yaml
import pandas as pd
from Bio import SeqIO


def load_config(config_file):
    with open(config_file, "r") as handle:
        return yaml.safe_load(handle)


def normalize_id(raw_id):
    """
    Normalize FASTA or table identifiers to core protein IDs.
    """

    if raw_id is None:
        return None

    x = str(raw_id).strip().split()[0]

    if x.startswith("transcript:"):
        x = x.replace("transcript:", "", 1)

    match = re.search(r"_cds_([A-Za-z]{2,}\d+(?:\.\d+)?)_", x)
    if match:
        return match.group(1)

    match = re.search(r"cds[-_:]([A-Za-z]{2,}\d+(?:\.\d+)?)", x)
    if match:
        return match.group(1)

    if "|" in x:
        for part in reversed(x.split("|")):
            match = re.search(r"([A-Za-z]{2,}\d+(?:\.\d+)?)", part)
            if match:
                return match.group(1)

    match = re.search(r"([A-Za-z]{2,}\d+(?:\.\d+)?)", x)
    if match:
        return match.group(1)

    return x


def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")


def inspect_fasta_ids(fasta_file, label, n=5):
    print(f"\n{label} ID normalization examples:")
    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        print(f"  {record.id} -> {normalize_id(record.id)}")
        if i + 1 >= n:
            break


def load_orthogroup_table(config):
    orthogroup_file = config["input"]["orthogroup_file"]
    check_file_exists(orthogroup_file)

    columns = config["columns"]

    df = pd.read_csv(orthogroup_file)

    required_cols = [
        columns["orthogroup_col"],
        columns["species1_protein_col"],
        columns["species2_protein_col"],
    ]

    missing_cols = [x for x in required_cols if x not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required columns in orthogroup table: {missing_cols}\n"
            f"Available columns: {list(df.columns)}"
        )

    df = df.rename(
        columns={
            columns["orthogroup_col"]: "Orthogroup",
            columns["species1_protein_col"]: "Species1_Protein",
            columns["species2_protein_col"]: "Species2_Protein",
        }
    )

    df = df.dropna(subset=["Orthogroup", "Species1_Protein", "Species2_Protein"])
    df["Species1_Protein_Normalized"] = df["Species1_Protein"].apply(normalize_id)
    df["Species2_Protein_Normalized"] = df["Species2_Protein"].apply(normalize_id)

    print(f"Loaded {len(df)} one-to-one orthologous gene pairs.")
    print("\nOrtholog table preview:")
    print(
        df[
            [
                "Orthogroup",
                "Species1_Protein",
                "Species1_Protein_Normalized",
                "Species2_Protein",
                "Species2_Protein_Normalized",
            ]
        ].head()
    )

    return df


def create_id_mapping(df, config):
    species1_code = config["species"]["species1_code"]
    species2_code = config["species"]["species2_code"]
    out_dir = config["output"]["work_dir"]

    os.makedirs(out_dir, exist_ok=True)

    mapping_rows = []

    for _, row in df.iterrows():
        og = row["Orthogroup"]

        mapping_rows.append(
            {
                "Orthogroup": og,
                "Original_ID": row["Species1_Protein"],
                "Normalized_ID": row["Species1_Protein_Normalized"],
                "New_ID": f"{species1_code}_{og}",
                "Species": species1_code,
            }
        )

        mapping_rows.append(
            {
                "Orthogroup": og,
                "Original_ID": row["Species2_Protein"],
                "Normalized_ID": row["Species2_Protein_Normalized"],
                "New_ID": f"{species2_code}_{og}",
                "Species": species2_code,
            }
        )

    mapping_df = pd.DataFrame(mapping_rows)

    mapping_file = os.path.join(out_dir, "id_mapping.tsv")
    mapping_df.to_csv(mapping_file, sep="\t", index=False)

    pairs_df = df[
        [
            "Orthogroup",
            "Species1_Protein_Normalized",
            "Species2_Protein_Normalized",
        ]
    ].copy()

    pairs_df.columns = [
        "Orthogroup",
        f"{species1_code}_Original",
        f"{species2_code}_Original",
    ]

    pairs_file = os.path.join(out_dir, "ortholog_pairs_clean.tsv")
    pairs_df.to_csv(pairs_file, sep="\t", index=False)

    print(f"\nID mapping table saved: {mapping_file}")
    print(f"Clean ortholog pair table saved: {pairs_file}")

    return mapping_df


def build_mapping_dict(mapping_df, species_code):
    sub = mapping_df[mapping_df["Species"] == species_code].copy()
    return dict(zip(sub["Normalized_ID"], sub["New_ID"]))


def rename_fasta_sequences(input_file, output_file, id_mapping, label, log_dir):
    print(f"\nRenaming {label}: {input_file}")

    total_count = 0
    renamed_count = 0
    unmatched = []

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            total_count += 1

            original_id = record.id
            normalized_id = normalize_id(original_id)

            if normalized_id in id_mapping:
                new_id = id_mapping[normalized_id]
                record.id = new_id
                record.name = new_id
                record.description = ""
                renamed_count += 1
            else:
                unmatched.append((original_id, normalized_id))

            SeqIO.write(record, out_handle, "fasta")

    print(f"  Total sequences: {total_count}")
    print(f"  Renamed sequences: {renamed_count}")
    print(f"  Unmatched sequences: {len(unmatched)}")
    print(f"  Output: {output_file}")

    if unmatched:
        log_file = os.path.join(log_dir, f"unmatched_ids_{label.replace(' ', '_')}.tsv")
        pd.DataFrame(unmatched, columns=["Original_ID", "Normalized_ID"]).to_csv(
            log_file, sep="\t", index=False
        )
        print(f"  Unmatched ID list saved: {log_file}")

    return renamed_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    species1_code = config["species"]["species1_code"]
    species2_code = config["species"]["species2_code"]

    renamed_dir = config["output"]["renamed_dir"]
    log_dir = config["output"]["logs_dir"]

    required_files = [
        config["input"]["orthogroup_file"],
        config["input"]["species1_cds"],
        config["input"]["species2_cds"],
        config["input"]["species1_protein"],
        config["input"]["species2_protein"],
    ]

    print("=" * 80)
    print("Step 1: Prepare ortholog data and rename FASTA sequences")
    print("=" * 80)

    for file_path in required_files:
        check_file_exists(file_path)
        print(f"Found: {file_path}")

    inspect_fasta_ids(config["input"]["species1_cds"], f"{species1_code} CDS")
    inspect_fasta_ids(config["input"]["species2_cds"], f"{species2_code} CDS")
    inspect_fasta_ids(config["input"]["species1_protein"], f"{species1_code} protein")
    inspect_fasta_ids(config["input"]["species2_protein"], f"{species2_code} protein")

    ortho_df = load_orthogroup_table(config)
    mapping_df = create_id_mapping(ortho_df, config)

    species1_map = build_mapping_dict(mapping_df, species1_code)
    species2_map = build_mapping_dict(mapping_df, species2_code)

    outputs = {
        "species1_cds": os.path.join(renamed_dir, f"{species1_code}_cds_renamed.fa"),
        "species2_cds": os.path.join(renamed_dir, f"{species2_code}_cds_renamed.fa"),
        "species1_protein": os.path.join(renamed_dir, f"{species1_code}_prot_renamed.faa"),
        "species2_protein": os.path.join(renamed_dir, f"{species2_code}_prot_renamed.faa"),
    }

    total_renamed = 0

    total_renamed += rename_fasta_sequences(
        config["input"]["species1_cds"],
        outputs["species1_cds"],
        species1_map,
        f"{species1_code}_CDS",
        log_dir,
    )

    total_renamed += rename_fasta_sequences(
        config["input"]["species2_cds"],
        outputs["species2_cds"],
        species2_map,
        f"{species2_code}_CDS",
        log_dir,
    )

    total_renamed += rename_fasta_sequences(
        config["input"]["species1_protein"],
        outputs["species1_protein"],
        species1_map,
        f"{species1_code}_protein",
        log_dir,
    )

    total_renamed += rename_fasta_sequences(
        config["input"]["species2_protein"],
        outputs["species2_protein"],
        species2_map,
        f"{species2_code}_protein",
        log_dir,
    )

    print("\nStep 1 completed successfully.")
    print(f"Ortholog pairs processed: {len(ortho_df)}")
    print(f"Total renamed sequences: {total_renamed}")


if __name__ == "__main__":
    main()
