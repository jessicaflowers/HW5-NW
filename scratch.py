#!/usr/bin/env python3
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


def main():
    #  small sequences you can compute by hand
    seqA = "GA"
    seqB = "GA"

    # Use BLOSUM62 and a linear gap penalty of -1 
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -1 )

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
