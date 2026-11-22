
from src.sequence_analysis import (
    read_fasta_file,
    calculate_sequence_statistics
)

from src.reference_database import load_reference_database
from src.gene_detection import detect_resistance_genes
from src.report_generation import (
    generate_report,
    display_report,
    export_report
)


def main():

    print("\nANTIBIOTIC RESISTANCE GENE DETECTION TOOL")
    print("----------------------------------------")

    file_path = input("\nEnter the path of your FASTA file: ")

    try:
        records = read_fasta_file(file_path)

        reference_genes = load_reference_database()

        all_results = []

        for record in records:

            print("\nSequence ID:", record["id"])

            statistics = calculate_sequence_statistics(
                record["sequence"]
            )

            print("Sequence length:", statistics["length"])
            print("GC percentage:", statistics["gc_percentage"], "%")

            results = detect_resistance_genes(
                record["sequence"],
                reference_genes
            )

            all_results.extend(results)

        report = generate_report(all_results)

        display_report(report)

        output_path = export_report(report)

        print("\nReport saved at:", output_path)

    except (ValueError, FileNotFoundError) as error:
        print("\nError:", error)


if __name__ == "__main__":
    main()
