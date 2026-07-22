#!/usr/bin/env python3

"""
check_protein_format.py

Check protein FASTA format and generate basic sequence statistics.
"""

import argparse
import os
from Bio import SeqIO


def check_protein_format(protein_file, log_file):
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")

    records = list(SeqIO.parse(protein_file, "fasta"))

    with open(log_file, "w") as log:
        log.write(f"Input file: {os.path.basename(protein_file)}\n")
        log.write("=" * 60 + "\n")

        if len(records) == 0:
            log.write("ERROR: No sequences were found.\n")
            return False

        log.write(f"Number of sequences: {len(records)}\n\n")

        log.write("Example sequence IDs:\n")
        for record in records[:5]:
            log.write(f"{record.id}\n")
        log.write("\n")

        lengths = [len(record.seq) for record in records]

        log.write("Sequence length statistics:\n")
        log.write(f"Mean length: {sum(lengths) / len(lengths):.2f}\n")
        log.write(f"Minimum length: {min(lengths)}\n")
        log.write(f"Maximum length: {max(lengths)}\n\n")

        invalid_chars = set()
        invalid_seq_count = 0

        for record in records:
            seq_chars = set(str(record.seq).upper())
            bad_chars = seq_chars - standard_aa
            if bad_chars:
                invalid_chars.update(bad_chars)
                invalid_seq_count += 1

        log.write("Amino acid composition check:\n")
        if invalid_chars:
            log.write(
                f"WARNING: {invalid_seq_count} sequences contain non-standard amino acid characters: "
                f"{', '.join(sorted(invalid_chars))}\n"
            )
        else:
            log.write("All sequences contain only 20 standard amino acids.\n")

        start_met_count = sum(
            1 for record in records
            if str(record.seq).upper().startswith("M")
        )

        log.write(
            f"\nSequences starting with methionine M: "
            f"{start_met_count} ({start_met_count / len(records) * 100:.2f}%)\n"
        )

        stop_count = sum(1 for record in records if "*" in str(record.seq))

        log.write(
            f"Sequences containing stop symbol '*': "
            f"{stop_count} ({stop_count / len(records) * 100:.2f}%)\n"
        )

    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input protein FASTA file")
    parser.add_argument("--log", required=True, help="Output log file")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    success = check_protein_format(args.input, args.log)

    if success:
        print(f"FASTA format check completed: {args.log}")
    else:
        raise RuntimeError(f"FASTA format check failed: {args.log}")


if __name__ == "__main__":
    main()
