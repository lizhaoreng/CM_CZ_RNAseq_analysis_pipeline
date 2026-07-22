#!/usr/bin/env python3

"""
clean_protein_fasta.py

Clean protein FASTA sequences before SignalP, TMHMM, and EffectorP prediction.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def clean_fasta(input_fasta, output_fasta):
    allowed = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ")
    records_out = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper()
        seq = seq.replace("*", "")
        seq = seq.replace(" ", "")
        seq = seq.replace("\n", "")
        seq = seq.replace("\r", "")
        seq = "".join([aa for aa in seq if aa in allowed])

        if len(seq) > 0:
            record.seq = Seq(seq)
            record.description = record.id
            records_out.append(record)

    SeqIO.write(records_out, output_fasta, "fasta")
    print(f"Cleaned sequences written: {len(records_out)}")
    print(f"Output file: {output_fasta}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input protein FASTA file")
    parser.add_argument("--output", required=True, help="Output cleaned FASTA file")
    args = parser.parse_args()

    clean_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
