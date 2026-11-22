
import streamlit as st
import pandas as pd

from src.sequence_analysis import (
    read_fasta_file,
    calculate_sequence_statistics
)

from src.reference_database import load_reference_database
from src.gene_detection import detect_resistance_genes
from src.report_generation import generate_report


st.set_page_config(
    page_title="Antibiotic Resistance Gene Detection",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Antibiotic Resistance Gene Detection Tool")

st.write(
    "Upload a bacterial DNA FASTA file to screen for "
    "known antibiotic resistance gene sequences."
)

st.info(
    "This tool performs sequence similarity screening. "
    "A detected sequence match does not confirm clinical resistance."
)

st.sidebar.header("Detection Settings")

identity_threshold = st.sidebar.slider(
    "Minimum percentage identity",
    min_value=50,
    max_value=100,
    value=90
)

coverage_threshold = st.sidebar.slider(
    "Minimum reference coverage",
    min_value=50,
    max_value=100,
    value=80
)

uploaded_file = st.file_uploader(
    "Upload a FASTA file",
    type=["fasta", "fa", "fna"]
)

if uploaded_file is not None:

    try:

        records = read_fasta_file(uploaded_file)

        reference_genes = load_reference_database()

        all_results = []
        statistics_data = []

        for record in records:

            statistics = calculate_sequence_statistics(
                record["sequence"]
            )

            statistics_data.append({
                "Sequence ID": record["id"],
                "Length (bp)": statistics["length"],
                "GC Content (%)": statistics["gc_percentage"]
            })

            results = detect_resistance_genes(
                record["sequence"],
                reference_genes,
                identity_threshold,
                coverage_threshold
            )

            for result in results:
                result["sequence_id"] = record["id"]

            all_results.extend(results)

        st.subheader("Sequence Statistics")

        statistics_df = pd.DataFrame(statistics_data)

        st.dataframe(
            statistics_df,
            use_container_width=True,
            hide_index=True
        )

        report = generate_report(all_results)

        if not report.empty:

            st.subheader("Detected Resistance Genes")

            st.dataframe(
                report,
                use_container_width=True,
                hide_index=True
            )

            st.success(
                f"{len(report)} qualifying gene match(es) detected."
            )

            csv_data = report.to_csv(index=False).encode("utf-8")

            st.download_button(
                label="Download Results as CSV",
                data=csv_data,
                file_name="detection_results.csv",
                mime="text/csv"
            )

        else:

            st.warning("No resistance genes detected.")

    except (ValueError, FileNotFoundError) as error:

        st.error(str(error))
