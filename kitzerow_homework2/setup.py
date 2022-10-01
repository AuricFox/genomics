# Samuel Kitzerow, kitze012
# Homework 1, Setting up matrices
# NOTE: not for main program, just for setting up lookup tables

import sys
import csv

# Matrix used to map codon inputs to index locations in the lists codon_count and codons
matrix = [
    [[0,  1,  2,  3], [4,  5,  6,  7], [8,  9,  10, 11], [12, 13, 14, 15]],     # Matrix A
    [[16, 17, 18, 19], [20, 21, 22, 23], [24, 25, 26, 27], [28, 29, 30, 31]],   # Matrix C
    [[32, 33, 34, 35], [36, 37, 38, 39], [40, 41, 42, 43], [44, 45, 46, 47]],   # Matrix G
    [[48, 49, 50, 51], [52, 53, 54, 55], [56, 57, 58, 59], [60, 61, 62, 63]]    # Matrix T
]
# Codon and amino acid pairs
data = [
    ["AAA",	"Lys"], ["AAC",	"Asn"], ["AAG",	"Lys"], ["AAT",	"Asn"], ["ACA",	"Thr"], ["ACC",	"Thr"], ["ACG",	"Thr"], ["ACT",	"Thr"],
    ["AGA",	"Arg"], ["AGC",	"Ser"], ["AGG",	"Arg"], ["AGT",	"Ser"], ["ATA",	"Ile"], ["ATC",	"Ile"], ["ATG",	"Met"], ["ATT",	"Ile"],
    ["CAA",	"Gln"], ["CAC",	"His"], ["CAG",	"Gln"], ["CAT",	"His"], ["CCA",	"Pro"], ["CCC",	"Pro"], ["CCG",	"Pro"], ["CCT",	"Pro"],
    ["CGA",	"Arg"], ["CGC",	"Arg"], ["CGG",	"Arg"], ["CGT",	"Arg"], ["CTA",	"Leu"], ["CTC",	"Leu"], ["CTG",	"Leu"], ["CTT",	"Leu"],
    ["GAA",	"Glu"], ["GAC",	"Asp"], ["GAG",	"Glu"], ["GAT",	"Asp"], ["GCA",	"Ala"],	["GCC",	"Ala"],	["GCG",	"Ala"],	["GCT",	"Ala"],	
    ["GGA",	"Gly"], ["GGC",	"Gly"], ["GGG",	"Gly"], ["GGT",	"Gly"], ["GTA",	"Val"], ["GTC",	"Val"], ["GTG",	"Val"], ["GTT",	"Val"],
    ["TAA",	"Stp"], ["TAC",	"Tyr"], ["TAG",	"Stp"], ["TAT",	"Tyr"], ["TCA",	"Ser"], ["TCC",	"Ser"], ["TCG",	"Ser"], ["TCT",	"Ser"],
    ["TGA",	"Stp"], ["TGC",	"Cys"], ["TGG",	"Trp"], ["TGT",	"Cys"], ["TTA",	"Leu"], ["TTC",	"Phe"], ["TTG",	"Leu"], ["TTT",	"Phe"]
]
# List of amino acids
amino = [  
    "Leu", "Val", "Thr", "Ala", "Ser", "Gly", "Lys", "Asn", "Ile", "Asp",
    "Phe", "Glu", "Tyr", "Pro", "Gln", "Arg", "Cys", "Met", "His", "Trp", "Stp"
]

# Matches the base character with an index in the matrix or returns -1
def char_to_index(base):
    if base == 'A':
        return 0
    elif base == 'C':
        return 1
    elif base == 'G':
        return 2
    elif base == 'T':
        return 3
    else:
        return -1

# ========================================================================================================================================
# Setting up lookup table that maps codons to amino acids
def moresetup():

    amatrix = []
    
    for i in range(64):     # Populating amatrix with zeros
        amatrix.append(0)

    for x in data:
        codon = x[0]
        m = char_to_index(codon[0])
        row = char_to_index(codon[1])
        column = char_to_index(codon[2])
        index = matrix[m][column][row]     # Finding the codon index in

        amatrix[index] = amino.index(x[1])

    print(amatrix)

# ========================================================================================================================================
# Merging csv files with corresponding values into one file
def merge_files(file1, file2, file3):
    f1 = open(file1, 'r')       # Separate genes
    f2 = open(file2, 'r')       # Whole genome

    data = []

    for (line1, line2) in zip(f1, f2):          # Iterates thru both files
        line1 = (line1.strip()).split(',')      # Removing newlines and splitting string at commas
        line2 = (line2.strip()).split(',')

        line1.append(line2[1])                  # Appending the count from the second file to a list from the first file
        data.append(line1)

    with open(file3, 'w', newline='') as file:
        cfile = csv.writer(file)
        cfile.writerows(data)                           # Wrights codon or amino acid data to csv
    

# ========================================================================================================================================
if __name__ == "__main__":

    # merge_files('SARS-CoV-2_separate_genes.csv', 'SARS-CoV-2_whole_genome.csv', 'combined_codons.csv')
    # merge_files('separate_amino_acids.csv', 'whole_amino_acids.csv', 'combined_amino_acids.csv')

    if(len(sys.argv) == 2 and sys.argv[1] == "setup"):      # Setting up amino acid matrix
        moresetup()

    elif(len(sys.argv) == 5 and sys.argv[1] == "merge"):    # Merging two csv files into one
        merge_files(sys.argv[2], sys.argv[3], sys.argv[4])

    
