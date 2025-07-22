from Bio import SeqIO
import sys

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(str(seq)))

def check_reverse_complements(forward_file, reverse_file):
    """
    Check if sequences in reverse_file are reverse complements of sequences in forward_file
    and in the same order.
    """
    # Read sequences from both files
    forward_seqs = list(SeqIO.parse(forward_file, "fasta"))
    reverse_seqs = list(SeqIO.parse(reverse_file, "fasta"))
    
    # Check if files have same number of sequences
    if len(forward_seqs) != len(reverse_seqs):
        print(f"❌ Files have different numbers of sequences:")
        print(f"   {forward_file}: {len(forward_seqs)} sequences")
        print(f"   {reverse_file}: {len(reverse_seqs)} sequences")
        return False
    
    all_match = True
    for i, (fwd, rev) in enumerate(zip(forward_seqs, reverse_seqs)):
        # Calculate expected reverse complement
        expected_rev = reverse_complement(fwd.seq)
        actual_rev = str(rev.seq)
        
        # Check if sequences match
        if expected_rev != actual_rev:
            print(f"\n❌ Mismatch at sequence {i+1}:")
            print(f"Forward ID: {fwd.id}")
            print(f"Reverse ID: {rev.id}")
            print(f"Expected reverse complement: {expected_rev[:50]}... (showing first 50 bases)")
            print(f"Actual reverse sequence:    {actual_rev[:50]}... (showing first 50 bases)")
            all_match = False
        
        # Check if IDs match (optional, might have different naming conventions)
        if fwd.id != rev.id:
            print(f"\nℹ️ Different sequence IDs at position {i+1}:")
            print(f"Forward ID: {fwd.id}")
            print(f"Reverse ID: {rev.id}")
    
    if all_match:
        print(f"\n✅ All sequences match their reverse complements!")
        print(f"Total sequences checked: {len(forward_seqs)}")
    
    return all_match

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_reverse_complements.py forward.fasta reverse.fasta")
        sys.exit(1)
        
    forward_file = sys.argv[1]
    reverse_file = sys.argv[2]
    check_reverse_complements(forward_file, reverse_file) 