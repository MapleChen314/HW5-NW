# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        self._align_matrix=np.zeros((len(self._seqA)+1,len(self._seqB)+1))
        self._gapA_matrix=np.tile(-np.inf,(len(self._seqA)+1,len(self._seqB)+1))
        self._gapB_matrix=np.tile(-np.inf,(len(self._seqA)+1,len(self._seqB)+1))
        self._back = np.tile("None",(len(self._seqA)+1,len(self._seqB)+1))
        #self._back[0,:]=np.tile(1,len(self._seqB)+1)
        self._back_A = np.tile("None",(len(self._seqA)+1,len(self._seqB)+1))
        self._back_B = np.tile("None",(len(self._seqA)+1,len(self._seqB)+1))
        for i in range(len(self._seqA)+1):
            self._gapA_matrix[i,0]=self.gap_open+self.gap_extend*i
            self._back_A[i,0]="eGA"
            self._align_matrix[i,0]=-np.inf
        for j in range(len(self._seqB)+1):
            self._gapB_matrix[0,j]=self.gap_open+self.gap_extend*j
            self._back_B[0,j]="eGB"
            self._align_matrix[0,j]=-np.inf
        self._align_matrix[0,0]=0
      
        # TODO: Implement global alignment here
        for i in range(1,len(self._seqA)+1):
            for j in range(1,len(self._seqB)+1):
                align_reward=self.sub_dict[(self._seqA[i-1],self._seqB[j-1])] #Sequences are 0-indexed and matrices are 1-indexed because of the extra column I have
                self._align_matrix[i,j]=max(
                    self._gapA_matrix[i-1,j-1]+align_reward, #Match from previous A gap
                    self._gapB_matrix[i-1,j-1]+align_reward, #Match from previous B gap
                    self._align_matrix[i-1,j-1]+align_reward #Match ; go diagonal
                )
                self._gapA_matrix[i,j]=max(
                    self._align_matrix[i-1,j]+self.gap_open+self.gap_extend, #Open gap in A; traverse vertically
                    self._gapA_matrix[i-1,j]+self.gap_extend #Extend gap
                )
                self._gapB_matrix[i,j]=max(
                    self._align_matrix[i,j-1]+self.gap_open+self.gap_extend, #open gap in B; traverse horizontally
                    self._gapB_matrix[i,j-1]+self.gap_extend #extend gap
                )
                back_idx=np.argmax(np.array((
                    self._gapA_matrix[i-1,j-1]+align_reward, 
                    self._gapB_matrix[i-1,j-1]+align_reward,
                    self._align_matrix[i-1,j-1]+align_reward 
                )))
                self._back[i,j]=["GA","GB","Mtch"][back_idx]
                back_A_idx=np.argmax(np.array((
                    self._align_matrix[i-1,j]+self.gap_open+self.gap_extend, #Open gap; go to _back and go vertical
                    self._gapA_matrix[i-1,j]+self.gap_extend #Extend gap; go to backA and go vertical
                )))
                self._back_A[i,j]=["oGA","eGA"][back_A_idx]
                back_B_idx=np.argmax(np.array((
                    self._align_matrix[i-1,j]+self.gap_open+self.gap_extend, #Open gap; go to _back and go horizontal
                    self._gapB_matrix[i-1,j]+self.gap_extend #Extend gap; go to backA and go horizontal
                )))
                self._back_B[i,j]=["oGB","eGB"][back_B_idx]
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        #Get the top scoring-- should this just be the bottom right?
        #[top_score_idxx,top_score_idxy]=np.unravel_index(np.argmax(self._align_matrix, axis=None), self._align_matrix.shape)
        self.alignment_score=max(
            self._align_matrix[-1,-1],
            self._gapA_matrix[-1,-1],
            self._gapB_matrix[-1,-1]
        )
        current_idx=[len(self._seqA),len(self._seqB)]
        current_matrix_idx=np.argmax(np.array((
            self._align_matrix[-1,-1],
            self._gapA_matrix[-1,-1],
            self._gapB_matrix[-1,-1])))
        matrices=[self._back,self._back_A,self._back_B]
        while (current_idx[0]>0 or current_idx[1]>0):
            pointer=matrices[current_matrix_idx][current_idx[0],current_idx[1]]
            if current_matrix_idx==0:
                if pointer=="GA": #Match from previous A Gap
                    change_idx=[-1,-1]
                    new_matrix_idx=1
                    self.seqA_align=self._seqA[current_idx[0]-1]+self.seqA_align
                    self.seqB_align=self._seqB[current_idx[1]-1]+self.seqB_align
                elif pointer=="GB": #Match from previous B Gap
                    change_idx=[-1,-1]
                    new_matrix_idx=2
                    self.seqA_align=self._seqA[current_idx[0]-1]+self.seqA_align
                    self.seqB_align=self._seqB[current_idx[1]-1]+self.seqB_align
                elif pointer=="Mtch": #Match from previous match
                    change_idx=[-1,-1]
                    new_matrix_idx=0
                    self.seqA_align=self._seqA[current_idx[0]-1]+self.seqA_align
                    self.seqB_align=self._seqB[current_idx[1]-1]+self.seqB_align
            elif current_matrix_idx==1: #Reference gapA
                change_idx=[-1,0]
                self.seqA_align=self._seqA[current_idx[0]-1]+self.seqA_align
                self.seqB_align="-"+self.seqB_align
                if pointer=="eGA":
                    new_matrix_idx=1
                else:
                    new_matrix_idx=0
            elif current_matrix_idx==2: #Reference GapB
                change_idx=[0,-1]
                self.seqA_align="-"+self.seqA_align
                self.seqB_align=self._seqB[current_idx[1]]+self.seqB_align
                if pointer=="eGA":
                    new_matrix_idx=2
                else:
                    new_matrix_idx=0
            current_idx=[x+y for x,y in zip(current_idx,change_idx)]
            current_matrix_idx=new_matrix_idx

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
