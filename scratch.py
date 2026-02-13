#!/usr/bin/env python3
"""
scratch.py

Tiny script to exercise the Needleman-Wunsch implementation on
very small sequences so the DP matrix and alignment can be checked by hand.

This script uses the linear-gap implementation in `align/align.py`.
"""

from align import NeedlemanWunsch
import numpy as np


def print_matrix(M: np.ndarray, seqA: str, seqB: str) -> None:
    """Pretty-print DP matrix M with sequence labels.

    M is expected to be shaped (len(seqA)+1, len(seqB)+1).
    """
    cols = M.shape[1]
    # Header
    header = [' '] + ['-'] + list(seqB)
    print('    ' + '  '.join(header))
    for i in range(M.shape[0]):
        row_label = '-' if i == 0 else seqA[i - 1]
        row_vals = [f"{M[i, j]:6.1f}" for j in range(cols)]
        print(f"{row_label} " + ' '.join(row_vals))


def main():
    # Very small sequences you can compute by hand
    seqA = "GA"
    seqB = "GA"

    # Use BLOSUM62 and a linear gap penalty of -1 (gap_extend ignored by linear implementation)
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -1 )#, -1)

    score, alnA, alnB = nw.align(seqA, seqB)

    # print("Sequence A:", seqA)
    # print("Sequence B:", seqB)
    # print()
    # print("DP matrix (M):")
    # print_matrix(nw._align_matrix, seqA, seqB)
    # print()
    # print("Alignment:")
    # print(alnA)
    # print(alnB)
    print("Score:", score)


if __name__ == '__main__':
    main()
