
import pytest

from src.sequence_analysis import (
    validate_dna_sequence,
    calculate_sequence_statistics,
    get_reverse_complement,
    read_fasta_file
)


def test_validate_dna_sequence():
    assert validate_dna_sequence("atgc") == "ATGC"


def test_invalid_dna_sequence():
    with pytest.raises(ValueError):
        validate_dna_sequence("ATGX")


def test_empty_sequence():
    with pytest.raises(ValueError):
        validate_dna_sequence("")


def test_sequence_statistics():
    result = calculate_sequence_statistics("ATGCATGC")

    assert result["length"] == 8
    assert result["A"] == 2
    assert result["T"] == 2
    assert result["G"] == 2
    assert result["C"] == 2
    assert result["gc_percentage"] == 50.0


def test_ambiguous_bases():
    result = calculate_sequence_statistics("ATGCN")

    assert result["length"] == 5
    assert result["ambiguous_bases"] == 1
    assert result["gc_percentage"] == 50.0


def test_reverse_complement():
    assert get_reverse_complement("ATGC") == "GCAT"


def test_read_fasta_file(tmp_path):
    fasta_file = tmp_path / "sample.fasta"

    fasta_file.write_text(
        ">sequence1\nATGCATGC\n"
        ">sequence2\nGGGGCCCC\n"
    )

    records = read_fasta_file(fasta_file)

    assert len(records) == 2
    assert records[0]["id"] == "sequence1"
    assert records[1]["sequence"] == "GGGGCCCC"


def test_missing_fasta_file():
    with pytest.raises(FileNotFoundError):
        read_fasta_file("missing_file.fasta")
