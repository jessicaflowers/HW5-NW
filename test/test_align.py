import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    # Use BLOSUM62 and a linear gap penalty of -10
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10)
    score, a, b = nw.align(seq1, seq2)

    # expected alignment matrix 
    expected_M = np.array([
        [0.0, -10.0, -20.0, -30.0],
        [-10.0, 5.0, -5.0, -15.0],
        [-20.0, -5.0, 4.0, -6.0],
        [-30.0, -15.0, 0.0, 5.0],
        [-40.0, -25.0, -10.0, 5.0],
    ])

    # back traced matrix expected values 
    expected_back = np.array([
        [0, 2, 2, 2],
        [1, 0, 2, 2],
        [1, 1, 0, 2],
        [1, 1, 0, 0],
        [1, 1, 1, 0],
    ], dtype=int)

    # Compare matrices
    assert nw._align_matrix.shape == expected_M.shape
    np.testing.assert_allclose(nw._align_matrix, expected_M, atol=1e-6)
    assert nw._back.shape == expected_back.shape
    np.testing.assert_array_equal(nw._back, expected_back)
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    # a gap penalty of -4 with BLOSUM62 should give a score of 18
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -4)
    score, a, b = nw.align(seq3, seq4)
    assert score == 18.0
    # expected alignments from README
    assert a == "MAVHQLIRRP"
    assert b == "M---QLIRHP"
    



