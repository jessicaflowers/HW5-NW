#!/usr/bin/env python3
from align import NeedlemanWunsch, read_fasta
from align import NeedlemanWunsch
import numpy as np


def print_matrix(M: np.ndarray, seqA: str, seqB: str) -> None:
    cols = M.shape[1]
    # Header
    header = [' '] + ['-'] + list(seqB)
    print('    ' + '  '.join(header))
    for i in range(M.shape[0]):
        row_label = '-' if i == 0 else seqA[i - 1]
        row_vals = [f"{M[i, j]:6.1f}" for j in range(cols)]
        print(f"{row_label} " + ' '.join(row_vals))



# #  small sequences you can compute by hand
# seqA = "GA"
# seqB = "GA"

# # Use BLOSUM62 and a linear gap penalty of -1 
# nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -1 )

# score, alnA, alnB = nw.align(seqA, seqB)

# # print("Sequence A:", seqA)
# # print("Sequence B:", seqB)
# # print()
# # print("DP matrix (M):")
# # print_matrix(nw._align_matrix, seqA, seqB)
# # print()
# # print("Alignment:")
# # print(alnA)
# # print(alnB)
# print("Score:", score)

seq1, _ = read_fasta("./data/test_seq1.fa")
seq2, _ = read_fasta("./data/test_seq2.fa")


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
print(nw._align_matrix)
# print(np.allclose(nw._align_matrix, expected_M, atol=1e-6))
print(nw._back)

# Compare matrices
# assert nw._align_matrix.shape == expected_M.shape
# np.testing.assert_allclose(nw._align_matrix, expected_M, atol=1e-6)
# assert nw._back.shape == expected_back.shape
# np.testing.assert_array_equal(nw._back, expected_back)
