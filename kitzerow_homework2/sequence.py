# Samuel Kitzerow, kitze012
# Homework 1, Counting Codons

# Codon types: A C G T
# Mapping codons to index values to an array/list used for counting
# 4 designator types with 3 total designators in a string (4^3 = 64 possible codon combinations)
# I manually coded in the codon combinations since the list was manageable otherwise I would 
# have used a tree or graph to generate the list of codons

# Matrix A:
# [
#     [AAA, ACA, AGA, ATA],   =>  [0,  1,  2,  3],
#     [AAC, ACC, AGC, ATC],   =>  [4,  5,  6,  7],
#     [AAG, ACG, AGG, ATG],   =>  [8,  9,  10, 11],
#     [AAT, ACT, AGT, ATT]    =>  [12, 13, 14, 15]
# ]

# Matrix C:
# [
#     [CAA, CCA, CGA, CTA],   =>  [16, 17, 18, 19]
#     [CAC, CCC, CGC, CTC],   =>  [20, 21, 22, 23]
#     [CAG, CCG, CGG, CTG],   =>  [24, 25, 26, 27]
#     [CAT, CCT, CGT, CTT]    =>  [28, 29, 30, 31]
# ]

# Matrix G:
# [
#     [GAA, GCA, GGA, GTA],   =>  [32, 33, 34, 35]
#     [GAC, GCC, GGC, GTC],   =>  [36, 37, 38, 39]
#     [GAG, GCG, GGG, GTG],   =>  [40, 41, 42, 43]
#     [GAT, GCT, GGT, GTT]    =>  [44, 45, 46, 47]
# ]

# Matrix T:
# [
#     [TAA, TCA, TGA, TTA],   =>  [48, 49, 50, 51]
#     [TAC, TCC, TGC, TTC],   =>  [52, 53, 54, 55]
#     [TAG, TCG, TGG, TTG],   =>  [56, 57, 58, 59]
#     [TAT, TCT, TGT, TTT]    =>  [60, 61, 62, 63]
# ]

class sequence:
    def __init__(self):
        self.codons = [ # List of all 64 possible DNA codon combinations in their corresponding index locations for codon_count
            "AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT",
            "CAA", "CCA", "CGA", "CTA", "CAC", "CCC", "CGC", "CTC", "CAG", "CCG", "CGG", "CTG", "CAT", "CCT", "CGT", "CTT",
            "GAA", "GCA", "GGA", "GTA", "GAC", "GCC", "GGC", "GTC", "GAG", "GCG", "GGG", "GTG", "GAT", "GCT", "GGT", "GTT",
            "TAA", "TCA", "TGA", "TTA", "TAC", "TCC", "TGC", "TTC", "TAG", "TCG", "TGG", "TTG", "TAT", "TCT", "TGT", "TTT"
        ]
        self.amino = [  # List of amino acids
            "Leu", "Val", "Thr", "Ala", "Ser", "Gly", "Lys", "Asn", "Ile", "Asp",
            "Phe", "Glu", "Tyr", "Pro", "Gln", "Arg", "Cys", "Met", "His", "Trp", "Stp"
        ]
        self.codon_count = [0] * 64   # Tracks the number of codons entered
        # Matrix used to map codon inputs to index locations in the lists codon_count and codons
        self.matrix = [
            [[0,  1,  2,  3], [4,  5,  6,  7], [8,  9,  10, 11], [12, 13, 14, 15]],     # Matrix A
            [[16, 17, 18, 19], [20, 21, 22, 23], [24, 25, 26, 27], [28, 29, 30, 31]],   # Matrix C
            [[32, 33, 34, 35], [36, 37, 38, 39], [40, 41, 42, 43], [44, 45, 46, 47]],   # Matrix G
            [[48, 49, 50, 51], [52, 53, 54, 55], [56, 57, 58, 59], [60, 61, 62, 63]]    # Matrix T
        ]
        # Amino acid lookup table used to map codons to amino acids (self.codons => self.amino)
        self.amatrix = [
            6, 2, 15, 8, 7, 2, 4, 8, 6, 2, 15, 17, 7, 2, 4, 8, 14, 13, 15, 0, 18, 13, 15, 0, 14, 13, 15, 0, 18, 13, 15, 0,
            11, 3, 5, 1, 9, 3, 5, 1, 11, 3, 5, 1, 9, 3, 5, 1, 20, 4, 20, 0, 12, 4, 16, 10, 20, 4, 19, 0, 12, 4, 16, 10
        ]

        self.sequence = ""          # Character string of bases
        self.header = ""            ## Character string for header info

        
    # ----------------------------------------------------------------------------------------------------------------------
    # Assigns the base character with an matrix/row/column value for reference in the matrix or returns -1 if invalid
    def char_to_index(self, base):
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

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the codon string into individual characters, maps them to the matching codon, and incriments the count of the codon
    def add_codon(self, codon):

        if (len(codon) != 3):                               # String can only have 3 bases
            print("ERROR: NOT A CODON!")
            return
        
        matrix = self.char_to_index(codon[0])
        row = self.char_to_index(codon[1])
        column = self.char_to_index(codon[2])

        if (matrix == -1 or row == -1 or column == -1):     # Incorrect char entered: A, C, G, T only
            print("ERROR: INVALID INPUT!")
            return

        index = self.matrix[matrix][column][row]            # Retreives mapped index in the matrix 
        # print("Index: ", index)
        self.codon_count[index] += 1                        # Incriments corresponding codon count

    # ----------------------------------------------------------------------------------------------------------------------
    # Adds codon or sequence strings to self.sequence
    def add_header(self, head):
        self.header = head

    # ----------------------------------------------------------------------------------------------------------------------
    # Adds codon or sequence strings to self.sequence
    def add_to_sequence(self, seq):
        self.sequence += seq

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the genome sequence string into individual codon sub-strings and adds them to the count
    def add_to_count(self, seq):
        end = len(seq) - len(seq) % 3   # Ignore any trailing characters

        for i in range(0, end, 3):      # Iterate thru each 3 chars in the string
            codon = seq[i:(i+3)]
            self.add_codon(codon)       # Add codon to count

    # ----------------------------------------------------------------------------------------------------------------------
    # Sums the totals of all codons and returns a list of codon and count pairs
    def get_codon_count(self):
        data = []
        
        for i in range(64):
            subdata = [self.codons[i], self.codon_count[i]]
            data.append(subdata)
    
        return data

    # ----------------------------------------------------------------------------------------------------------------------
    # Sums total of all amino acids and returns a list of amino acid and count pairs
    def get_amino_count(self):
        amino_count = []

        for i in range(len(self.amino)):            # Populating amino_count with zeros
            amino_count.append([self.amino[i], 0])

        for x in range(len(self.codon_count)):      # Adding up the amino acids
            count = self.codon_count[x]             # Number of codons counted
            index = self.amatrix[x]                 # Index of where to add the count

            amino_count[index][1] += count          # Updating amino acid count

        return amino_count

    # ----------------------------------------------------------------------------------------------------------------------
    # Testing if elements are being populated properly
    def debugging(self):
        for i in range(64):
            print(self.codons[i], self.codon_count[i], sep=" : ")


if __name__ == "__main__":

    test = sequence()
    codons = [ # list of all 64 possible DNA codon combinations
        "AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT",
        "CAA", "CCA", "CGA", "CTA", "CAC", "CCC", "CGC", "CTC", "CAG", "CCG", "CGG", "CTG", "CAT", "CCT", "CGT", "CTT",
        "GAA", "GCA", "GGA", "GTA", "GAC", "GCC", "GGC", "GTC", "GAG", "GCG", "GGG", "GTG", "GAT", "GCT", "GGT", "GTT",
        "TAA", "TCA", "TGA", "TTA", "TAC", "TCC", "TGC", "TTC", "TAG", "TCG", "TGG", "TTG", "TAT", "TCT", "TGT", "TTT"
    ]

    for x in codons:
        test.add_codon(x)
    
    # test.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    test.debugging()
    # print(test.get_amino_count())
    