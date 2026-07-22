#!/usr/bin/env python3

"""
extract_secretome_candidates.py

Extract protein sequences based on SignalP6 or TMHMM output.

Modes:
  signalp: extract SignalP-positive proteins from SignalP6 prediction_results.txt
  tmhmm:   extract proteins with no predicted transmembrane helices from TMHMM output
"""

import argparse
import re
from Bio import SeqIO


def extract_signalp_positive(signalp_txt):
    positive_ids = set()

    with open(signalp_txt) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = re.split(r"\s+", line)

            if len(parts) < 2:
                continue

            protein_id = parts[0]
            fields = parts[1:]

            if any(x == "SP" or x.startswith("SP(") for x in fields):
                positive_ids.add(protein_id)

    return positive_ids


def extract_tmhmm_no_tm(tmhmm_txt):
    no_tm_ids = set()

    with open(tmhmm_txt) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            # Short TMHMM format:
            # ID len=xxx ExpAA=... First60=... PredHel=0 Topology=o
            match_short = re.search(r"PredHel=(\d+)", line)

            if match_short:
                protein_id = line.split()[0]
                predhel = int(match_short.group(1))

                if predhel == 0:
                    no_tm_ids.add(protein_id)

                continue

            # Long TMHMM format:
            # # ID Number of predicted TMHs:  0
            match_long = re.match(
                r"^#\s+(.+?)\s+Number of predicted TMHs:\s+(\d+)",
                line
            )

            if match_long:
                protein_id = match_long.group(1)
                predhel = int(match_long.group(2))

                if predhel == 0:
                    no_tm_ids.add(protein_id)

    return no_tm_ids


def write_selected_sequences(input_fasta, selected_ids, output_fasta):
    records = []
    total = 0

    for record in SeqIO.parse(input_fasta, "fasta"):
        total += 1

        if record.id in selected_ids:
            records.append(record)

    SeqIO.write(records, output_fasta, "fasta")

    print(f"Input FASTA sequences: {total}")
    print(f"Selected IDs: {len(selected_ids)}")
    print(f"Sequences written: {len(records)}")
    print(f"Output FASTA: {output_fasta}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", required=True, help="SignalP or TMHMM result file")
    parser.add_argument("--fasta", required=True, help="Input protein FASTA")
    parser.add_argument("--output", required=True, help="Output selected FASTA")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["signalp", "tmhmm"],
        help="Extraction mode"
    )

    args = parser.parse_args()

    if args.mode == "signalp":
        selected_ids = extract_signalp_positive(args.result)
    elif args.mode == "tmhmm":
        selected_ids = extract_tmhmm_no_tm(args.result)
    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

    if len(selected_ids) == 0:
        print(f"WARNING: No selected IDs found for mode: {args.mode}")

    write_selected_sequences(args.fasta, selected_ids, args.output)


if __name__ == "__main__":
    main()
