from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from tabulate import tabulate

from Bio import SeqIO


def load_sequence(filename):
    """Load a sequence from a FASTA file and return it as a string"""
    record = next(SeqIO.parse(filename, "fasta"))
    return str(record.seq)


print("✅ All libraries working!")
def find_snps(ref, seq):
    snps = []
    for i, (r, s) in enumerate(zip(ref, seq)):
        if r != s and r in "ATGC" and s in "ATGC":
            snps.append((i, r, s))
    return snps

def calculate_ts_tv(snps):
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    ts = 0
    tv = 0

    for _, ref_base, alt_base in snps:
        if (ref_base, alt_base) in transitions:
            ts += 1
        elif ref_base != alt_base:  
            tv += 1

    ratio = ts / tv if tv > 0 else float("inf")
    return ts, tv, ratio

import matplotlib.pyplot as plt

def plot_snp_distribution(snps, title):
    positions = [pos for pos, _, _ in snps]

    plt.figure(figsize=(10, 3))
    plt.scatter(positions, [1]*len(positions), marker="|", color="red")
    plt.title(title)
    plt.xlabel("Position in reference sequence")
    plt.yticks([])  
    plt.show()

def build_frequency_matrix(snps):
    bases = ["A", "T", "C", "G"]
    matrix = pd.DataFrame(0.0, index=bases, columns=bases)

    
    for _, ref, alt in snps:
        matrix.loc[alt, ref] += 1

   
    total = matrix.values.sum()
    if total > 0:
        matrix = matrix / total

    return matrix




def format_ascii_table(matrix):
    bases = ["A", "T", "C", "G"]
    headers = ["Base after SNP \\ Original Base"] + bases
    table_data = []
    for row_base in bases:
        row = [row_base]
        for col_base in bases:
            if row_base == col_base:
                row.append("X")
            else:
                row.append(matrix.loc[row_base, col_base])
        table_data.append(row)
    return tabulate(table_data, headers=headers, tablefmt="grid")


ref = load_sequence("reference.fasta")
seq1 = load_sequence("sequence_1.fasta")
seq2 = load_sequence("sequence_2.fasta")

snps1 = find_snps(ref, seq1)
matrix1 = build_frequency_matrix(snps1)

print("=== SNP Frequency Matrix (ASCII style) ===")
print(format_ascii_table(matrix1))

snps2 = find_snps(ref, seq2)
matrix2 = build_frequency_matrix(snps2)

print("\n=== SNP Frequency Matrix for sequence_2.fasta ===")
print(format_ascii_table(matrix2))

plot_snp_distribution(snps1, "SNP Distribution: sequence_1 vs reference")
plot_snp_distribution(snps2, "SNP Distribution: sequence_2 vs reference")

ts1, tv1, ratio1 = calculate_ts_tv(snps1)
print(f"\nSample 1: Transitions={ts1}, Transversions={tv1}, Ts/Tv={ratio1:.2f}")

ts2, tv2, ratio2 = calculate_ts_tv(snps2)
print(f"Sample 2: Transitions={ts2}, Transversions={tv2}, Ts/Tv={ratio2:.2f}")


