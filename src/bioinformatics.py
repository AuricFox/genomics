# GUI terminal
# Runs main program that calls other scripts

import sys
import os
import sequence as sq
import alignment as al
import de_bruijn as db
import variance as vr
import utils

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp")

# ==============================================================================================================
"""
Extracts the header and sequence data from the input file and counts the totals in the sequence
Parameter(s):
    * file(str): path to the file containing the sequence data
Output:
    * res(dict): dictionary response containing the totals of codons, amino acids, and/or kmers
"""
# ==============================================================================================================
def get_sequence_data(file:str, codon=None, amino=None, kmer=None, k:int=3):
    data = utils.get_data(os.path.join(path, file))
    res = {}

    if codon != None:
        seq = sq.Codon(data[0], "Codons")
        res['codons'] = seq.codon
    if amino != None:
        seq = sq.Amino_Acid(data[0], "Amino Acids")
        res['amino acids'] = seq.amino_acid
    if kmer != None:
        seq = sq.Kmer(data[0], "K-mers", k)
        res['kmers'] = seq.kmers

    return res

# ==============================================================================================================
"""
Extracts the header and sequence data from the input file and counts the totals in the sequence
Parameter(s):
    * file1 (str): path to the file containing the first sequence data
    * file2 (str): path to the file containing the second sequence data
    * gap_pen (int): penalty for gaps in the alignment
    * match_point (int): point(s) added to the score for matches in the alignment
    * match_pen (int): penalty for mismatches in the alignment
    * ignore (bool): condition to ignore start and end gap penalties for local alignment
Output:
    * res(dict): dictionary response containing the totals of codons, amino acids, and/or kmers
"""
# ==============================================================================================================
def get_alignment_data(file1:str, file2:str, gap_pen:int=-2, match_point:int=1, match_pen:int=-1, ignore:bool=False):
    seq1 = utils.get_data(os.path.join(path, file1))
    seq2 = utils.get_data(os.path.join(path, file2))

    data = al.Alignment(
        seq=seq1[0][0], 
        ref=seq2[0][0], 
        gap_pen=gap_pen, 
        match_point=match_point, 
        match_pen=match_pen, 
        ignore=ignore
    )

    res = data.results

    return res
# ==============================================================================================================
def main():
    data = get_alignment_data(file1='pfizer_mrna.fna', file2='sars_spike_protein.fna',  ignore=True)
    return


if __name__ == "__main__":
    main()
