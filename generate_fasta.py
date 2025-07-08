import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def _resolve_output_path(output_filename: str, same_dir_as_script: bool) -> str:
    if same_dir_as_script:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_path = os.path.join(script_dir, output_filename)
    else:
        output_path = output_filename
    # Make sure the folder exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    return output_path


def generate_fasta(seq: str, output_filename: str, seq_id: str, description: str = "", same_dir_as_script=True):
    output_path = _resolve_output_path(output_filename, same_dir_as_script)
    record = SeqRecord(Seq(seq), id=seq_id, description=description)
    SeqIO.write(record, output_path, "fasta")
    print(f"FASTA written to {output_path}")


def generate_mult_fasta(seqs, output_filename, same_dir_as_script=True):
    """
    seqs: list of tuples (sequence_str, seq_id, [description])
    output_filename: name of the FASTA file to write
    """
    output_path = _resolve_output_path(output_filename, same_dir_as_script)

    records = []
    for entry in seqs:
        if len(entry) == 3:
            seq_str, seq_id, description = entry
        elif len(entry) == 2:
            seq_str, seq_id = entry
            description = ""
        else:
            raise ValueError("Each item in `seqs` must be a tuple of (sequence_str, seq_id, [description])")

        record = SeqRecord(Seq(seq_str), id=seq_id, description=description)
        records.append(record)

    SeqIO.write(records, output_path, "fasta")
    print(f"Wrote {len(records)} sequences to {output_path}")

    """
    Seq formatting for generate_mult_fasta(): 
    
    seqs = [
        ("ATGCGT", "seq1", "example sequence 1"),
        ("ATGAAA", "seq2", "example sequence 2"),
        ("ATGTTT", "seq3", "example sequence 3"),
    ]
    """



