#!/usr/bin/env python3

"""
02_extract_ortholog_pairs.py

Extract paired protein and CDS sequences for one-to-one orthologous groups.
"""

import argparse
import os
import sys
import yaml
import pandas as pd
from Bio import SeqIO


def load_config(config_file):
    with open(config_file, "r") as handle:
        return yaml.safe_load(handle)


def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Required file not found: {file_path}")


def load_sequences(fasta_file, label):
    print(f"Loading {label}: {fasta_file}")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    print(f"  Loaded {len(seq_dict)} sequences.")
    return seq_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    species1_code = config["species"]["species1_code"]
    species2_code = config["species"]["species2_code"]

    work_dir = config["output"]["work_dir"]
    renamed_dir = config["output"]["renamed_dir"]
    pair_prot_dir = config["output"]["pair_protein_dir"]
    pair_cds_dir = config["output"]["pair_cds_dir"]
    log_dir = config["output"]["logs_dir"]

    mapping_file = os.path.join(work_dir, "id_mapping.tsv")

    species1_cds = os.path.join(renamed_dir, f"{species1_code}_cds_renamed.fa")
    species2_cds = os.path.join(renamed_dir, f"{species2_code}_cds_renamed.fa")
    species1_protein = os.path.join(renamed_dir, f"{species1_code}_prot_renamed.faa")
    species2_protein = os.path.join(renamed_dir, f"{species2_code}_prot_renamed.faa")

    print("=" * 80)
    print("Step 2: Extract one-to-one orthologous sequence pairs")
    print("=" * 80)

    for file_path in [mapping_file, species1_cds, species2_cds, species1_protein, species2_protein]:
        check_file_exists(file_path)
        print(f"Found: {file_path}")

    os.makedirs(pair_prot_dir, exist_ok=True)
    os.makedirs(pair_cds_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    mapping_df = pd.read_csv(mapping_file, sep="\t")
    orthogroups = sorted(mapping_df["Orthogroup"].unique())

    print(f"\nNumber of orthogroups: {len(orthogroups)}")

    sp1_prot = load_sequences(species1_protein, f"{species1_code} protein")
    sp2_prot = load_sequences(species2_protein, f"{species2_code} protein")
    sp1_cds = load_sequences(species1_cds, f"{species1_code} CDS")
    sp2_cds = load_sequences(species2_cds, f"{species2_code} CDS")

    successful_pairs = 0
    missing_rows = []

    for index, og in enumerate(orthogroups, start=1):
        sp1_id = f"{species1_code}_{og}"
        sp2_id = f"{species2_code}_{og}"

        missing_items = []

        if sp1_id not in sp1_prot:
            missing_items.append(f"{sp1_id}:protein")
        if sp2_id not in sp2_prot:
            missing_items.append(f"{sp2_id}:protein")
        if sp1_id not in sp1_cds:
            missing_items.append(f"{sp1_id}:CDS")
        if sp2_id not in sp2_cds:
            missing_items.append(f"{sp2_id}:CDS")

        if missing_items:
            missing_rows.append(
                {"Orthogroup": og, "Missing": ";".join(missing_items)}
            )
            continue

        prot_file = os.path.join(pair_prot_dir, f"{og}.faa")
        cds_file = os.path.join(pair_cds_dir, f"{og}.cds")

        with open(prot_file, "w") as handle:
            handle.write(f">{sp1_id}\n{str(sp1_prot[sp1_id].seq)}\n")
            handle.write(f">{sp2_id}\n{str(sp2_prot[sp2_id].seq)}\n")

        with open(cds_file, "w") as handle:
            handle.write(f">{sp1_id}\n{str(sp1_cds[sp1_id].seq)}\n")
            handle.write(f">{sp2_id}\n{str(sp2_cds[sp2_id].seq)}\n")

        successful_pairs += 1

        if successful_pairs % 500 == 0:
            print(f"  Extracted {successful_pairs} sequence pairs.")

    if missing_rows:
        missing_file = os.path.join(log_dir, "missing_sequences.tsv")
        pd.DataFrame(missing_rows).to_csv(missing_file, sep="\t", index=False)
        print(f"\nMissing sequence details saved: {missing_file}")

    print("\nStep 2 completed.")
    print(f"Successful pairs: {successful_pairs}")
    print(f"Missing pairs: {len(missing_rows)}")
    print(f"Protein pair directory: {pair_prot_dir}")
    print(f"CDS pair directory: {pair_cds_dir}")

    if successful_pairs == 0:
        print("ERROR: No valid sequence pairs were extracted.")
        sys.exit(1)


if __name__ == "__main__":
    main()
