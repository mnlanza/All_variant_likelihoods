#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os.path

def read_fasta(fasta_file):
    """Read a FASTA file and return a dictionary of id -> sequence"""
    sequences = {}
    # Expand the user path before parsing
    expanded_path = os.path.expanduser(fasta_file)
    for record in SeqIO.parse(expanded_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def compare_fastas(file1, file2):
    """Compare two FASTA files and report differences"""
    print(f"Comparing {file1} with {file2}")
    
    # Read both files
    seqs1 = read_fasta(file1)
    seqs2 = read_fasta(file2)
    
    # Compare sequence IDs
    ids1 = set(seqs1.keys())
    ids2 = set(seqs2.keys())
    
    print(f"\nFile 1 has {len(ids1)} sequences")
    print(f"File 2 has {len(ids2)} sequences")
    
    # Check for unique IDs in each file
    unique_to_1 = ids1 - ids2
    unique_to_2 = ids2 - ids1
    
    if unique_to_1:
        print(f"\nSequences only in {file1}:")
        for seq_id in unique_to_1:
            print(f"  {seq_id}")
    
    if unique_to_2:
        print(f"\nSequences only in {file2}:")
        for seq_id in unique_to_2:
            print(f"  {seq_id}")
    
    # Compare common sequences
    common_ids = ids1 & ids2
    print(f"\nComparing {len(common_ids)} common sequences...")
    
    differences = []
    for seq_id in common_ids:
        if seqs1[seq_id] != seqs2[seq_id]:
            differences.append(seq_id)
            len1 = len(seqs1[seq_id])
            len2 = len(seqs2[seq_id])
            print(f"\nDifference in {seq_id}:")
            print(f"  Length in file1: {len1}")
            print(f"  Length in file2: {len2}")
            if len1 == len2:
                # Find position of first difference
                pos = next(i for i in range(len1) if seqs1[seq_id][i] != seqs2[seq_id][i])
                print(f"  First difference at position {pos + 1}")
                context = 10  # Show 10 bases before and after difference
                start = max(0, pos - context)
                end = min(len1, pos + context + 1)
                print(f"  File1: ...{seqs1[seq_id][start:end]}...")
                print(f"  File2: ...{seqs2[seq_id][start:end]}...")
    
    if not differences:
        print("\nAll common sequences are identical!")
    else:
        print(f"\nFound differences in {len(differences)} sequences")
        print("Sequences with differences:")
        for seq_id in differences:
            print(f"  {seq_id}")

if __name__ == "__main__":
    file1 = "~/Documents/margin_human/input/human_contigs_src.fasta"
    file2 = "~/Documents/margin_human/input/human_contigs_src_test.fasta"
    compare_fastas(file1, file2) 