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
def get_sequence_data(file:str, codon:str=None, amino:str=None, kmer:str=None, k:int=3):
    '''
    Extracts the header and sequence data from the input file and counts the totals in the sequence
    
    Parameter(s):
        file (str): path to the file containing the sequence data
        codon (str, optional): get the codon count if not None
        amino (str, optional): get the amino acid count if not None
        kmer (str, optional): get the k-mer count if not None
        k (int, optional): size of the k-mer
    
    Output:
        response (dict): dictionary response containing the totals of codons, amino acids, and/or kmers
    '''

    data = utils.get_data(os.path.join(path, file))
    response = {}

    if codon != None:
        seq = sq.Codon(sequences=data[0], header="Codons")
        response['codons'] = seq.codon
    if amino != None:
        seq = sq.Amino_Acid(sequences=data[0], header="Amino Acids")
        response['amino acids'] = seq.amino_acid
    if kmer != None:
        seq = sq.Kmer(sequences=data[0], header="K-mers", k=k)
        response['kmers'] = seq.kmers

    return response

# ==============================================================================================================
def get_alignment_data(file1:str, file2:str, gap_pen:int=-2, match_point:int=1, match_pen:int=-1, ignore:bool=False):
    '''
    Extracts the header and sequence data from the input file and counts the totals in the sequence
    
    Parameter(s):
        file1 (str): path to the file containing the first sequence data
        file2 (str): path to the file containing the second sequence data
        gap_pen (int, optional): penalty for gaps in the alignment
        match_point (int, optional): point(s) added to the score for matches in the alignment
        match_pen (int, optional): penalty for mismatches in the alignment
        ignore (bool, optional): condition to ignore start and end gap penalties for local alignment
    
    Output:
        response (dict): dictionary response containing the totals of codons, amino acids, and/or kmers
    '''

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

    response = data.results

    return response

# ==============================================================================================================
def get_variance_data(file:str, plot:bool, smooth:bool, vRegion:bool, spacing:int=30, numV:int=6,):
    '''
    Calculates the variance between the sequences and plots the data
    
    Parameter(s):
        file (str): path to the file containing the sequence data
        plot (bool): save the raw/smooth plot to a file if true
        smooth (bool): smooth the data if true, else plot the raw data
        vRegion (bool): save the plot with variance regions to a file if true
        spacing (int, optional): specifies the minimum number of bases needed for a v region (peak)
        numV (int, optional): specifies the number of desired v regions (peaks)
    
    Output(s):
    '''

    seq = utils.get_data(os.path.join(path, file))
    data = vr.Variance(sequences=seq[0], header=seq[1])

    data.save_plot(
        plot=plot, 
        smooth=smooth, 
        vRegion=vRegion, 
        spacing=spacing, 
        numV=numV
    )

    return True

# ==============================================================================================================
def main():
    data = get_variance_data(file='sequences.fna', plot=True, smooth=True, vRegion=True)
    return


if __name__ == "__main__":
    main()
