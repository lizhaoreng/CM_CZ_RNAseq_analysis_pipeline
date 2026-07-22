#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate pairwise p-distances from a concatenated FASTA alignment.

Usage:
    python 06_calculate_pairwise_distance.py \
        concatenated/Cercospora_ITS_TEF1_ACT_CAL_HIS3_concat.fasta \
        pairwise_distance_5gene.tsv

If output file is omitted, results are printed to stdout.
"""

import sys
from itertools import combinations


def read_fasta(path):
    seqs = {}
    name = None
    buf = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()

                name = line[1:].split()[0]
                buf = []

            else:
                buf.append(line)

        if name is not None:
            seqs[name] = "".join(buf).upper()

    return seqs


def pdist(seq1, seq2):
    comp = 0
    diff = 0

    missing = set("-?Nn")

    for x, y in zip(seq1, seq2):
        if x in missing or y in missing:
            continue

        comp += 1

        if x != y:
            diff += 1

    p = diff / comp if comp > 0 else float("nan")

    return diff, comp, p


def fmt_p(p):
    if p != p:
        return "nan"

    return f"{p:.6f}"


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(
            "Usage: python 06_calculate_pairwise_distance.py "
            "<concat_alignment.fasta> [output.tsv]\n"
        )
        sys.exit(1)

    fasta = sys.argv[1]
    outpath = sys.argv[2] if len(sys.argv) >= 3 else None

    seqs = read_fasta(fasta)

    targets = [
        "CM_Jilin_this_study",
        "CZ_Jilin_this_study",
        "C_zeae_maydis_SCOH1_5_genome",
        "C_zeina_CMW25467_genome",
        "C_zeae_maydis_2197",
        "C_zeae_maydis_CBS117763",
        "C_zeae_maydis_Czm_20_107",
        "C_zeina_CPC11995_extype",
        "C_zeina_CPC11998",
        "C_zeina_CMW25467",
    ]

    targets_present = [x for x in targets if x in seqs]
    missing_targets = [x for x in targets if x not in seqs]

    lines = []
    lines.append("taxon1\ttaxon2\tdiff\tcompared\tp_distance")

    for a, b in combinations(targets_present, 2):
        d, c, p = pdist(seqs[a], seqs[b])
        lines.append(f"{a}\t{b}\t{d}\t{c}\t{fmt_p(p)}")

    text = "\n".join(lines) + "\n"

    if outpath:
        with open(outpath, "w") as out:
            out.write(text)
    else:
        print(text, end="")

    sys.stderr.write(f"Loaded sequences: {len(seqs)}\n")
    sys.stderr.write(f"Target taxa present: {len(targets_present)}\n")

    if missing_targets:
        sys.stderr.write("Missing target taxa:\n")
        for x in missing_targets:
            sys.stderr.write(f"  {x}\n")

    key_pairs = [
        ("CM_Jilin_this_study", "C_zeae_maydis_SCOH1_5_genome"),
        ("CM_Jilin_this_study", "C_zeina_CMW25467_genome"),
        ("CZ_Jilin_this_study", "C_zeina_CMW25467_genome"),
        ("CZ_Jilin_this_study", "C_zeae_maydis_SCOH1_5_genome"),
        ("C_zeae_maydis_SCOH1_5_genome", "C_zeina_CMW25467_genome"),
    ]

    sys.stderr.write("\nKey comparisons:\n")
    sys.stderr.write("taxon1\ttaxon2\tdiff\tcompared\tp_distance\n")

    for a, b in key_pairs:
        if a in seqs and b in seqs:
            d, c, p = pdist(seqs[a], seqs[b])
            sys.stderr.write(f"{a}\t{b}\t{d}\t{c}\t{fmt_p(p)}\n")
        else:
            sys.stderr.write(f"{a}\t{b}\tNA\tNA\tNA\n")


if __name__ == "__main__":
    main()
