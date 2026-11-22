
import pandas as pd
from Bio import SeqIO
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

FASTA_FILE = BASE_DIR / "data" / "reference_genes.fasta"
INFO_FILE = BASE_DIR / "data" / "gene_information.csv"


def load_reference_database():

    if not FASTA_FILE.exists():
        raise FileNotFoundError("Reference FASTA file was not found.")

    if not INFO_FILE.exists():
        raise FileNotFoundError("Gene information file was not found.")

    sequences = list(SeqIO.parse(FASTA_FILE, "fasta"))
    information = pd.read_csv(INFO_FILE)

    if len(sequences) == 0:
        raise ValueError("No reference sequences found.")

    reference_genes = []

    for record in sequences:

        gene_info = information[
            information["gene_id"] == record.id
        ]

        if gene_info.empty:
            raise ValueError("Missing information for gene: " + record.id)

        gene = gene_info.iloc[0]

        reference_genes.append({
            "gene_id": gene["gene_id"],
            "gene_name": gene["gene_name"],
            "antibiotic_class": gene["antibiotic_class"],
            "mechanism": gene["mechanism"],
            "accession": gene["accession"],
            "source": gene["source"],
            "sequence": str(record.seq)
        })

    return reference_genes
