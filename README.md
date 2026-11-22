
# Antibiotic Resistance Gene Detection and Analysis Tool

## 1. Project Overview

The Antibiotic Resistance Gene Detection and Analysis Tool is a Python-based bioinformatics application designed to identify known antibiotic resistance genes in bacterial DNA sequences.

The tool compares uploaded DNA sequences against a small reference database using local sequence alignment. It also calculates sequence statistics, reports matching genes and exports detection results as CSV files.

## 2. Objectives

- Identify known antibiotic resistance gene sequences.
- Analyse DNA sequence length and GC content.
- Calculate percentage identity and reference coverage.
- Detect sequences in forward and reverse-complement orientations.
- Generate structured detection reports.
- Provide a simple web interface for FASTA file analysis.

## 3. Technologies Used

- Python
- Biopython
- Pandas
- Streamlit
- Pytest
- Visual Studio Code

## 4. Reference Gene Database

The application uses three reference resistance genes obtained from NCBI:

| Gene | Accession | Length |
|---|---|---:|
| blaTEM-1 | EF035599.1 | 861 bp |
| tetA | L06940.1 | 1218 bp |
| sulI | AY224185.1 | 840 bp |

The reference sequences were obtained from annotated GenBank records.

## 5. Methodology

1. Upload bacterial DNA sequences in FASTA format.
2. Validate the nucleotide sequences.
3. Calculate sequence length, nucleotide counts and GC content.
4. Load the reference resistance gene database.
5. Perform local alignment using Biopython PairwiseAligner.
6. Analyse forward and reverse-complement orientations.
7. Calculate percentage identity and reference coverage.
8. Report matches that satisfy the selected thresholds.
9. Export detection results as a CSV file.

## 6. Detection Parameters

Default screening thresholds:

- Minimum percentage identity: 90%
- Minimum reference coverage: 80%

These are prototype screening thresholds and are not universal clinical criteria.

## 7. Application Features

- FASTA file upload
- Multiple sequence processing
- DNA sequence validation
- GC content calculation
- Resistance gene screening
- Forward and reverse-complement detection
- Adjustable detection thresholds
- Results table with sequence IDs
- CSV report download

## 8. Testing

The project includes 12 automated tests using Pytest.

All 12 automated tests passed during project testing.

Additional manual tests include:

- Exact reference sequence detection
- Reverse-complement detection
- Unrelated sequence rejection
- Synthetic mutated sequence detection
- Positive FASTA sample
- Negative FASTA sample
- Multiple FASTA sequence processing
- CSV report generation

## 9. Project Limitations

- The reference database currently contains only three resistance genes.
- Detection results depend on the selected identity and coverage thresholds.
- Sequence similarity alone does not establish functional or clinical antibiotic resistance.
- The application is intended for educational and research screening, not clinical diagnosis.

## 10. Future Improvements

- Expand the reference database with additional resistance genes.
- Integrate more comprehensive resistance gene databases.
- Improve alignment and mutation analysis.
- Add more detailed visual reports.
- Evaluate the tool using larger and independently sourced bacterial datasets.

## 11. How to Run the Project

### Install the required libraries

```bash
pip install -r requirements.txt
```

### Run the web application

```bash
streamlit run app.py
```

Open the local URL displayed in the terminal to access the application.

### Run automated tests

```bash
python -m pytest -v
```

## 12. Author

Developed as an academic Computational Biology project.

## 13. GitHub Repository

[Antibiotic Resistance Gene Detection and Analysis Tool](https://github.com/sarayudevarakondamu-crypto/Antibiotic-Resistance-Gene-Detection)

