
from src.reference_database import load_reference_database
from src.gene_detection import detect_resistance_genes
from src.sequence_analysis import get_reverse_complement


def test_exact_sequence_detection():

    genes = load_reference_database()
    query = genes[0]["sequence"]

    results = detect_resistance_genes(query, genes)

    assert len(results) == 1
    assert results[0]["gene_name"] == "blaTEM-1"
    assert results[0]["identity"] == 100.0
    assert results[0]["coverage"] == 100.0


def test_reverse_complement_detection():

    genes = load_reference_database()
    query = get_reverse_complement(genes[0]["sequence"])

    results = detect_resistance_genes(query, genes)

    assert len(results) == 1
    assert results[0]["gene_name"] == "blaTEM-1"
    assert results[0]["strand"] == "reverse"


def test_unrelated_sequence():

    genes = load_reference_database()
    query = "ATGC" * 200

    results = detect_resistance_genes(query, genes)

    assert len(results) == 0


def test_mutated_sequence_detection():

    genes = load_reference_database()
    original = genes[0]["sequence"]

    changes = {
        "A": "C",
        "C": "A",
        "G": "T",
        "T": "G"
    }

    query = "".join(
        changes[base] if i % 20 == 0 else base
        for i, base in enumerate(original)
    )

    results = detect_resistance_genes(query, genes)

    assert len(results) == 1
    assert results[0]["gene_name"] == "blaTEM-1"
    assert results[0]["identity"] >= 90.0
    assert results[0]["coverage"] >= 80.0
