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
# Codon dictionary 
codon_dict = {
    "AAA":	{"amino_acid": "Lys", "count": 0}, "AAC":	{"amino_acid": "Asn", "count": 0}, "AAG":	{"amino_acid": "Lys", "count": 0}, 
    "AAT":	{"amino_acid": "Asn", "count": 0}, "ACA":	{"amino_acid": "Thr", "count": 0}, "ACC":	{"amino_acid": "Thr", "count": 0}, 
    "ACG":	{"amino_acid": "Thr", "count": 0}, "ACT":	{"amino_acid": "Thr", "count": 0}, "AGA":	{"amino_acid": "Arg", "count": 0}, 
    "AGC":	{"amino_acid": "Ser", "count": 0}, "AGG":	{"amino_acid": "Arg", "count": 0}, "AGT":	{"amino_acid": "Ser", "count": 0}, 
    "ATA":	{"amino_acid": "Ile", "count": 0}, "ATC":	{"amino_acid": "Ile", "count": 0}, "ATG":	{"amino_acid": "Met", "count": 0}, 
    "ATT":	{"amino_acid": "Ile", "count": 0}, "CAA":	{"amino_acid": "Gln", "count": 0}, "CAC":	{"amino_acid": "His", "count": 0}, 
    "CAG":	{"amino_acid": "Gln", "count": 0}, "CAT":	{"amino_acid": "His", "count": 0}, "CCA":	{"amino_acid": "Pro", "count": 0}, 
    "CCC":	{"amino_acid": "Pro", "count": 0}, "CCG":	{"amino_acid": "Pro", "count": 0}, "CCT":	{"amino_acid": "Pro", "count": 0},
    "CGA":	{"amino_acid": "Arg", "count": 0}, "CGC":	{"amino_acid": "Arg", "count": 0}, "CGG":	{"amino_acid": "Arg", "count": 0}, 
    "CGT":	{"amino_acid": "Arg", "count": 0}, "CTA":	{"amino_acid": "Leu", "count": 0}, "CTC":	{"amino_acid": "Leu", "count": 0}, 
    "CTG":	{"amino_acid": "Leu", "count": 0}, "CTT":	{"amino_acid": "Leu", "count": 0}, "GAA":	{"amino_acid": "Glu", "count": 0}, 
    "GAC":	{"amino_acid": "Asp", "count": 0}, "GAG":	{"amino_acid": "Glu", "count": 0}, "GAT":	{"amino_acid": "Asp", "count": 0}, 
    "GCA":	{"amino_acid": "Ala", "count": 0}, "GCC":	{"amino_acid": "Ala", "count": 0}, "GCG":	{"amino_acid": "Ala", "count": 0},	
    "GCT":	{"amino_acid": "Ala", "count": 0}, "GGA":	{"amino_acid": "Gly", "count": 0}, "GGC":	{"amino_acid": "Gly", "count": 0}, 
    "GGG":	{"amino_acid": "Gly", "count": 0}, "GGT":	{"amino_acid": "Gly", "count": 0}, "GTA":	{"amino_acid": "Val", "count": 0}, 
    "GTC":	{"amino_acid": "Val", "count": 0}, "GTG":	{"amino_acid": "Val", "count": 0}, "GTT":	{"amino_acid": "Val", "count": 0},
    "TAA":	{"amino_acid": "Stp", "count": 0}, "TAC":	{"amino_acid": "Tyr", "count": 0}, "TAG":	{"amino_acid": "Stp", "count": 0}, 
    "TAT":	{"amino_acid": "Tyr", "count": 0}, "TCA":	{"amino_acid": "Ser", "count": 0}, "TCC":	{"amino_acid": "Ser", "count": 0}, 
    "TCG":	{"amino_acid": "Ser", "count": 0}, "TCT":	{"amino_acid": "Ser", "count": 0}, "TGA":	{"amino_acid": "Stp", "count": 0}, 
    "TGC":	{"amino_acid": "Cys", "count": 0}, "TGG":	{"amino_acid": "Trp", "count": 0}, "TGT":	{"amino_acid": "Cys", "count": 0}, 
    "TTA":	{"amino_acid": "Leu", "count": 0}, "TTC":	{"amino_acid": "Phe", "count": 0}, "TTG":	{"amino_acid": "Leu", "count": 0}, 
    "TTT":	{"amino_acid": "Phe", "count": 0}
}
# Amino dictionary
amino_dict = { 
    "Leu": {"letter": "L", "count": 0}, "Val": {"letter": "V", "count": 0}, "Thr": {"letter": "T", "count": 0}, "Ala": {"letter": "A", "count": 0}, 
    "Ser": {"letter": "S", "count": 0}, "Gly": {"letter": "G", "count": 0}, "Lys": {"letter": "K", "count": 0}, "Asn": {"letter": "N", "count": 0}, 
    "Ile": {"letter": "I", "count": 0}, "Asp": {"letter": "D", "count": 0}, "Phe": {"letter": "F", "count": 0}, "Glu": {"letter": "E", "count": 0}, 
    "Tyr": {"letter": "Y", "count": 0}, "Pro": {"letter": "P", "count": 0}, "Gln": {"letter": "Q", "count": 0}, "Arg": {"letter": "R", "count": 0}, 
    "Cys": {"letter": "C", "count": 0}, "Met": {"letter": "M", "count": 0}, "His": {"letter": "H", "count": 0}, "Trp": {"letter": "W", "count": 0}, 
    "Stp": {"letter": "O", "count": 0}
}

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

    
