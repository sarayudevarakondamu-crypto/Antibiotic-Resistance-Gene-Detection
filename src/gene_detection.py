
from Bio.Align import PairwiseAligner

from src.sequence_analysis import (
    validate_dna_sequence,
    get_reverse_complement
)


def create_aligner():

    aligner = PairwiseAligner()

    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    return aligner


def calculate_alignment_statistics(alignment, reference_length):

    coordinates = alignment.coordinates

    matches = 0
    alignment_length = 0
    reference_bases_aligned = 0

    for i in range(len(coordinates[0]) - 1):

        ref_start = int(coordinates[0][i])
        ref_end = int(coordinates[0][i + 1])

        query_start = int(coordinates[1][i])
        query_end = int(coordinates[1][i + 1])

        ref_difference = ref_end - ref_start
        query_difference = query_end - query_start

        if ref_difference > 0 and query_difference > 0:

            ref_part = alignment.target[ref_start:ref_end]
            query_part = alignment.query[query_start:query_end]

            for j in range(ref_difference):
                if ref_part[j] == query_part[j]:
                    matches += 1

            alignment_length += ref_difference
            reference_bases_aligned += ref_difference

        elif ref_difference > 0:

            alignment_length += ref_difference
            reference_bases_aligned += ref_difference

        elif query_difference > 0:

            alignment_length += query_difference

    if alignment_length > 0:
        identity = (matches / alignment_length) * 100
    else:
        identity = 0.0

    if reference_length > 0:
        coverage = (reference_bases_aligned / reference_length) * 100
    else:
        coverage = 0.0

    return round(identity, 2), round(coverage, 2)


def detect_single_gene(query_sequence, reference_gene,
                       aligner, identity_threshold, coverage_threshold):

    reference_sequence = validate_dna_sequence(
        reference_gene["sequence"]
    )

    query_sequence = validate_dna_sequence(query_sequence)

    best_result = None

    for strand, sequence in [
        ("forward", query_sequence),
        ("reverse", get_reverse_complement(query_sequence))
    ]:

        alignments = aligner.align(reference_sequence, sequence)

        if len(alignments) == 0:
            continue

        alignment = alignments[0]

        identity, coverage = calculate_alignment_statistics(
            alignment,
            len(reference_sequence)
        )

        if identity < identity_threshold:
            continue

        if coverage < coverage_threshold:
            continue

        query_start = int(alignment.aligned[1][0][0])
        query_end = int(alignment.aligned[1][-1][1])

        if strand == "reverse":
            original_start = len(query_sequence) - query_end
            original_end = len(query_sequence) - query_start
        else:
            original_start = query_start
            original_end = query_end

        result = {
            "gene_id": reference_gene["gene_id"],
            "gene_name": reference_gene["gene_name"],
            "antibiotic_class": reference_gene["antibiotic_class"],
            "mechanism": reference_gene["mechanism"],
            "identity": identity,
            "coverage": coverage,
            "start_position": original_start + 1,
            "end_position": original_end,
            "strand": strand,
            "accession": reference_gene["accession"],
            "source": reference_gene["source"]
        }

        if best_result is None:
            best_result = result

        elif identity * coverage > best_result["identity"] * best_result["coverage"]:
            best_result = result

    return best_result


def detect_resistance_genes(query_sequence, reference_genes,
                            identity_threshold=90.0,
                            coverage_threshold=80.0):

    query_sequence = validate_dna_sequence(query_sequence)

    if not 0 <= identity_threshold <= 100:
        raise ValueError("Identity threshold must be between 0 and 100.")

    if not 0 <= coverage_threshold <= 100:
        raise ValueError("Coverage threshold must be between 0 and 100.")

    aligner = create_aligner()

    detected_genes = []

    for reference_gene in reference_genes:

        result = detect_single_gene(
            query_sequence,
            reference_gene,
            aligner,
            identity_threshold,
            coverage_threshold
        )

        if result is not None:
            detected_genes.append(result)

    return detected_genes
