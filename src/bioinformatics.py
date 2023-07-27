# GUI terminal
# Runs main program that calls other scripts

import sys
import os
import sequence as sq
import alignment as al
import de_bruijn as db
import utils

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp")

# ==============================================================================================================
# Extracts the header and sequence data from the input file and counts the totals in the sequence
# Input:
#   * file(str): path to the file containing the sequence data
# Returns:
#   * res(dict): dictionary response containing the totals of codons, amino acids, and/or kmers
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
# Extracts the header and sequence data from the input file and counts the totals in the sequence
# Input:
#   * file(str): path to the file containing the sequence data
# Returns:
#   * res(dict): dictionary response containing the totals of codons, amino acids, and/or kmers
# ==============================================================================================================
def get_alignment_data(file1:str, file2:str, option:str):
    seq1 = utils.get_data(os.path.join(path, file1))
    seq2 = utils.get_data(os.path.join(path, file2))
    res = {}

    

    return res
# ==============================================================================================================
def main():
    data = get_sequence_data('testing.fna', kmer='kmer', k=4)
    print(data)
    return


if __name__ == "__main__":
    main()
