from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")
    # use BLOSUM62 and linear gap penalty of -10 
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10) 

    def get_species_name(header: str) -> str:
        start = header.index('OS=') + 3 # the name is between OS ans OX in the data files
        end = header.index(' OX=', start)
        return header[start:end]

    species = [
        (get_species_name(gg_header), gg_seq),
        (get_species_name(mm_header), mm_seq),
        (get_species_name(br_header), br_seq),
        (get_species_name(tt_header), tt_seq),
    ]

    results = []
    for name, seq in species:
        score, alnA, alnB = nw.align(hs_seq, seq)
        results.append((name, float(score), alnA, alnB))

    # sort by score descending (most similar first)
    results.sort(key=lambda x: x[1], reverse=True)

    # Align all species to humans and print species in order of most similar to human BRD
    print("species ordered by similarity to human BRD2:")
    for i, (name, score, a, b) in enumerate(results, start=1):
        print(f"{i}. {name}: {score}")

    # Print all alignment scores between each species BRD2 and human BRD2
    print('\n')
    print("alignments and scores:")
    for name, score, a, b in results:
        print(f"{name} and human alignment score: {score}")
        # print(a)
        # print(b)
    

if __name__ == "__main__":
    main()
