"""Basic test examples."""

from mito_pipeline.pipeline import check_external_tools
from mito_pipeline.pipeline import get_softclip_lengths, get_max_indel_length


def test_check_external_tools_returns_list():
    # Regardless of whether minimap2/samtools/medaka are installed in this
    # environment, this function should still run and return a list.
    result = check_external_tools()
    assert isinstance(result, list)


def test_get_softclip_lengths_none_cigar():
    assert get_softclip_lengths(None) == 0


def test_get_softclip_lengths_sums_softclip_ops():
    # CIGAR op 4 = soft-clip
    cigartuples = [(4, 10), (0, 100), (4, 5)]
    assert get_softclip_lengths(cigartuples) == 15


def test_get_max_indel_length():
    # CIGAR op 1 = insertion, op 2 = deletion
    cigartuples = [(0, 100), (1, 3), (2, 7), (0, 50)]
    assert get_max_indel_length(cigartuples) == 7


def test_get_max_indel_length_no_indels():
    cigartuples = [(0, 100)]
    assert get_max_indel_length(cigartuples) == 0


# TODO: add integration tests for find_best_reference() and
# run_medaka_iterative() (needs small test FASTQ/FASTA fixtures, and an
# environment with minimap2/samtools/medaka installed)
