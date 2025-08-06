#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os.path
from typing import Dict, Tuple, List

def read_fasta(fasta_file: str) -> Dict[str, str]:
    """Read a FASTA file and return a dictionary of id -> sequence"""
    sequences = {}
    expanded_path = os.path.expanduser(fasta_file)
    for record in SeqIO.parse(expanded_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def find_sequence_difference(seq1: str, seq2: str) -> List[Tuple[int, str, str]]:
    """Find all positions where two sequences differ"""
    differences = []
    if len(seq1) != len(seq2):
        return differences
    
    for i, (base1, base2) in enumerate(zip(seq1, seq2)):
        if base1 != base2:
            differences.append((i + 1, base1, base2))  # 1-based position
    return differences

def compare_fastas_with_mutations(file1: str, file2: str) -> None:
    """
    Compare two FASTA files and check for mutations.
    Handles cases where one file might have double the sequences (original + mutated variants).
    """
    print(f"Comparing {file1} with {file2}")
    
    # Read both files
    seqs1 = read_fasta(file1)
    seqs2 = read_fasta(file2)
    
    # Get sequence counts
    count1 = len(seqs1)
    count2 = len(seqs2)
    
    print(f"\nFile 1 has {count1} sequences")
    print(f"File 2 has {count2} sequences")
    
    # Check if one file has double the sequences
    is_double = (count1 == 2 * count2) or (count2 == 2 * count1)
    
    if is_double:
        print("\nDetected that one file has double the sequences of the other!")
        # Determine which file has more sequences
        if count1 > count2:
            shorter_seqs = seqs2
            longer_seqs = seqs1
            shorter_file = file2
            longer_file = file1
        else:
            shorter_seqs = seqs1
            longer_seqs = seqs2
            shorter_file = file1
            longer_file = file2
            
        print("\nAnalyzing first sequence from shorter file against first two sequences in longer file...")
        
        # Get first sequence from shorter file
        first_short_id = list(shorter_seqs.keys())[0]
        first_short_seq = shorter_seqs[first_short_id]
        
        # Get first two sequences from longer file
        first_two_long_ids = list(longer_seqs.keys())[:2]
        
        # Compare sequences
        for long_id in first_two_long_ids:
            long_seq = longer_seqs[long_id]
            differences = find_sequence_difference(first_short_seq, long_seq)
            
            if differences:
                print(f"\nComparing {first_short_id} with {long_id}:")
                print(f"Found {len(differences)} difference(s):")
                for pos, base1, base2 in differences:
                    print(f"  Position {pos}: {base1} → {base2}")
                
                # Show context around first difference
                if differences:
                    pos = differences[0][0] - 1  # Convert to 0-based
                    context = 10
                    start = max(0, pos - context)
                    end = min(len(first_short_seq), pos + context + 1)
                    print(f"\nContext around first difference:")
                    print(f"  {first_short_id}: .......{first_short_seq[start:end]}...")
                    print(f"  {long_id}: ...{long_seq[start:end]}...")
            else:
                print(f"\n{first_short_id} and {long_id} are identical!")
    
    else:
        # Regular comparison for non-double cases
        ids1 = set(seqs1.keys())
        ids2 = set(seqs2.keys())
        
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
        
        common_ids = ids1 & ids2
        print(f"\nComparing {len(common_ids)} common sequences...")
        
        differences_found = False
        for seq_id in common_ids:
            differences = find_sequence_difference(seqs1[seq_id], seqs2[seq_id])
            if differences:
                differences_found = True
                print(f"\nDifferences in {seq_id}:")
                for pos, base1, base2 in differences:
                    print(f"  Position {pos}: {base1} → {base2}")
        
        if not differences_found:
            print("\nAll common sequences are identical!")

if __name__ == "__main__":
    # Example usage with your FASTA files
    file1 = "~/Documents/margin_human/input/human_contigs_src.fasta"
    file2 = "~/Documents/layer_human/input/mut_human_src_v1.fasta"
    compare_fastas_with_mutations(file1, file2) 