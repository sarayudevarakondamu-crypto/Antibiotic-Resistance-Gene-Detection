
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from io import StringIO


VALID_DNA_CHARACTERS = "ACGTNRYKMSWBDHV"


def validate_dna_sequence(sequence):

    sequence = str(sequence).upper()
    sequence = "".join(sequence.split())

    if len(sequence) == 0:
        raise ValueError("The DNA sequence is empty.")

    for base in sequence:
        if base not in VALID_DNA_CHARACTERS:
            raise ValueError("Invalid DNA character found: " + base)

    return sequence


def read_fasta_file(file_path):

    records = []

    if hasattr(file_path, "read"):

        content = file_path.getvalue().decode("utf-8")

        fasta_records = SeqIO.parse(
            StringIO(content),
            "fasta"
        )

    else:

        path = Path(file_path)

        if not path.exists():
            raise FileNotFoundError("FASTA file was not found.")

        if not path.is_file():
            raise ValueError("The given path is not a file.")

        if path.stat().st_size == 0:
            raise ValueError("The FASTA file is empty.")

        fasta_records = SeqIO.parse(str(path), "fasta")

    for record in fasta_records:

        sequence = validate_dna_sequence(record.seq)

        sequence_info = {
            "id": record.id,
            "description": record.description,
            "sequence": sequence
        }

        records.append(sequence_info)

    if len(records) == 0:
        raise ValueError("No valid FASTA sequences were found.")

    return records


def calculate_sequence_statistics(sequence):

    dna = validate_dna_sequence(sequence)

    a_count = dna.count("A")
    t_count = dna.count("T")
    g_count = dna.count("G")
    c_count = dna.count("C")

    total_bases = a_count + t_count + g_count + c_count
    gc_count = g_count + c_count

    if total_bases > 0:
        gc_percentage = (gc_count / total_bases) * 100
    else:
        gc_percentage = 0.0

    statistics = {
        "length": len(dna),
        "A": a_count,
        "T": t_count,
        "G": g_count,
        "C": c_count,
        "ambiguous_bases": len(dna) - total_bases,
        "gc_percentage": round(gc_percentage, 2)
    }

    return statistics


def get_reverse_complement(sequence):

    dna = validate_dna_sequence(sequence)

    reverse_complement = Seq(dna).reverse_complement()

    return str(reverse_complement)
