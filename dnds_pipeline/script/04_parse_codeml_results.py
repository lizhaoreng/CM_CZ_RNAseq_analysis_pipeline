#!/usr/bin/env python3

"""
04_parse_codeml_results.py

Parse pairwise codeml output files, merge alignment quality statistics,
apply quality-control filters, and generate clean dN/dS result tables.
"""

import argparse
import os
import re
import yaml
import pandas as pd
import numpy as np


def load_config(config_file):
    with open(config_file, "r") as handle:
        return yaml.safe_load(handle)


def parse_codeml_output(output_file):
    with open(output_file, "r") as handle:
        content = handle.read()

    patterns = {
        "omega": [
            r"omega\s*\(dN/dS\)\s*=\s*([0-9.+\-eE]+)",
            r"dN/dS\s*=\s*([0-9.+\-eE]+)",
            r"\bw\s*=\s*([0-9.+\-eE]+)",
            r"omega\s*=\s*([0-9.+\-eE]+)",
        ],
        "dN": [
            r"\bdN\s*=\s*([0-9.+\-eE]+)",
            r"\bdN\s*:\s*([0-9.+\-eE]+)",
        ],
        "dS": [
            r"\bdS\s*=\s*([0-9.+\-eE]+)",
            r"\bdS\s*:\s*([0-9.+\-eE]+)",
        ],
        "kappa": [
            r"kappa\s*\(ts/tv\)\s*=\s*([0-9.+\-eE]+)",
            r"kappa\s*=\s*([0-9.+\-eE]+)",
        ],
        "lnL": [
            r"lnL.*?=\s*([-0-9.+\-eE]+)",
        ],
        "tree_length": [
            r"tree length\s*=\s*([0-9.+\-eE]+)",
        ],
    }

    result = {}

    for key, regex_list in patterns.items():
        result[key] = np.nan
        for regex in regex_list:
            match = re.search(regex, content, re.IGNORECASE)
            if match:
                try:
                    value = float(match.group(1))
                    result[key] = value
                    break
                except ValueError:
                    continue

    return result


def add_original_ids(df, mapping_file):
    if not os.path.exists(mapping_file):
        print(f"Warning: ID mapping file not found: {mapping_file}")
        return df

    mapping_df = pd.read_csv(mapping_file, sep="\t")

    for species in mapping_df["Species"].unique():
        species_map = (
            mapping_df[mapping_df["Species"] == species]
            .set_index("Orthogroup")["Original_ID"]
            .to_dict()
        )
        df[f"{species}_Original_ID"] = df["Orthogroup"].map(species_map)

    return df


def apply_quality_control(df, qc_config):
    original_count = len(df)

    df["QC_Pass"] = True
    df["QC_Reason"] = "Pass"

    df["Has_Internal_Stop"] = df["Internal_Stops"].fillna(0) > 0

    filters = [
        (df["dS"].isna(), "dS_missing"),
        (df["dN"].isna(), "dN_missing"),
        (df["omega"].isna(), "omega_missing"),
        (df["Effective_Codons"] < qc_config["min_effective_codons"], "too_short_codons"),
        (df["Gap_Percentage"] > qc_config["max_gap_percentage"], "too_many_gaps"),
        (df["Sequence_Identity"] < qc_config["min_sequence_identity"], "low_sequence_identity"),
        (df["dS"] == 0, "dS_zero"),
        (df["dS"] < 0.005, "dS_too_low"),
        (df["dS"] > qc_config["max_dS"], "dS_too_high"),
        (df["omega"] <= 0, "omega_non_positive"),
        (df["omega"] > qc_config["max_omega"], "omega_extreme"),
        (df["Has_Internal_Stop"], "internal_stop_codon"),
        (df["dN"] > 5, "dN_extreme"),
    ]

    print("\nQuality-control filtering:")

    for condition, reason in filters:
        mask = condition & df["QC_Pass"]
        df.loc[mask, "QC_Pass"] = False
        df.loc[mask, "QC_Reason"] = reason

        if mask.sum() > 0:
            print(f"  {reason}: {mask.sum()} genes filtered")

    passed = df["QC_Pass"].sum()
    print(f"\nQC result: {passed}/{original_count} genes passed ({passed / original_count * 100:.1f}%)")

    return df


def classify_selection(omega):
    if pd.isna(omega):
        return "Unknown"
    if omega > 1:
        return "Accelerated (omega > 1)"
    if abs(omega - 1) < 0.05:
        return "Near neutral"
    return "Purifying (omega < 1)"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    work_dir = config["output"]["work_dir"]
    codeml_dir = config["output"]["codeml_out_dir"]

    length_stats_file = os.path.join(work_dir, "alignment_length_stats.csv")
    mapping_file = os.path.join(work_dir, "id_mapping.tsv")

    complete_out = os.path.join(work_dir, "dnds_enhanced_complete.csv")
    clean_out = os.path.join(work_dir, "dnds_enhanced_clean.csv")
    accelerated_out = os.path.join(work_dir, "accelerated_genes.csv")

    print("=" * 80)
    print("Step 4: Parse codeml results and apply quality control")
    print("=" * 80)

    if not os.path.exists(length_stats_file):
        raise FileNotFoundError(f"Missing alignment statistics file: {length_stats_file}")

    if not os.path.exists(codeml_dir):
        raise FileNotFoundError(f"codeml output directory not found: {codeml_dir}")

    length_stats = pd.read_csv(length_stats_file)

    results = []

    output_files = sorted([f for f in os.listdir(codeml_dir) if f.endswith(".out")])
    print(f"Found {len(output_files)} codeml output files.")

    for i, filename in enumerate(output_files, start=1):
        orthogroup = filename.replace(".out", "")
        file_path = os.path.join(codeml_dir, filename)

        parsed = parse_codeml_output(file_path)
        parsed["Orthogroup"] = orthogroup
        results.append(parsed)

        if i % 1000 == 0:
            print(f"Parsed {i}/{len(output_files)} files.")

    if not results:
        raise RuntimeError("No codeml results were parsed.")

    dnds_df = pd.DataFrame(results)
    dnds_df = dnds_df.merge(length_stats, on="Orthogroup", how="left")

    dnds_df = apply_quality_control(dnds_df, config["quality_filter"])
    dnds_df["Selection_Type"] = dnds_df["omega"].apply(classify_selection)

    dnds_df = add_original_ids(dnds_df, mapping_file)

    dnds_df.to_csv(complete_out, index=False)

    clean_df = dnds_df[dnds_df["QC_Pass"]].copy()
    clean_df.to_csv(clean_out, index=False)

    accelerated = clean_df[clean_df["omega"] > 1].sort_values("omega", ascending=False)
    accelerated.to_csv(accelerated_out, index=False)

    print("\nOutputs saved:")
    print(f"  Complete table: {complete_out}")
    print(f"  QC-passed table: {clean_out}")
    print(f"  Accelerated genes: {accelerated_out}")

    print("\nSummary of QC-passed omega values:")
    print(clean_df["omega"].describe())

    print("\nSelection category distribution:")
    print(clean_df["Selection_Type"].value_counts())


if __name__ == "__main__":
    main()
