
import pandas as pd
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def generate_report(detected_genes):

    columns = [
        "gene_id",
        "gene_name",
        "antibiotic_class",
        "mechanism",
        "identity",
        "coverage",
        "start_position",
        "end_position",
        "strand",
        "accession",
        "source"
    ]

    report = pd.DataFrame(detected_genes, columns=columns)

    return report


def export_report(report, filename="detection_results.csv"):

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    output_path = RESULTS_DIR / filename

    report.to_csv(output_path, index=False)

    return output_path


def display_report(report):

    if report.empty:
        print("No resistance genes detected.")
    else:
        print("\nResistance Gene Detection Results")
        print(report.to_string(index=False))
