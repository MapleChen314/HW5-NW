# Import NeedlemanWunsch class and read_fasta function
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

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    NW=NeedlemanWunsch(sub_matrix_file='./substitution_matrices/BLOSUM62.mat',gap_open=-10,gap_extend=-1)
    seqs=[gg_seq,mm_seq,br_seq,tt_seq]
    names=["Gallus gallus","Mus musculus","Balaeniceps rex"," Tursiops truncatus"]
    score_dict={}
    for i in range(len(seqs)):
        seq=seqs[i]
        name=names[i]
        NW.align(hs_seq,seq)
        score_dict[name]=NW.alignment_score

    ranked_dict=dict(sorted(score_dict.items(), key=lambda item: item[1],reverse=True))
    for name in ranked_dict.keys():
        print(name)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for score in ranked_dict.values():
        print(score)
    

if __name__ == "__main__":
    main()
