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

class map_codons:
    def __init__(self):
        self.codons = [ # list of all 64 possible DNA codon combinations in their corresponding index locations for codon_count
            "AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT",
            "CAA", "CCA", "CGA", "CTA", "CAC", "CCC", "CGC", "CTC", "CAG", "CCG", "CGG", "CTG", "CAT", "CCT", "CGT", "CTT",
            "GAA", "GCA", "GGA", "GTA", "GAC", "GCC", "GGC", "GTC", "GAG", "GCG", "GGG", "GTG", "GAT", "GCT", "GGT", "GTT",
            "TAA", "TCA", "TGA", "TTA", "TAC", "TCC", "TGC", "TTC", "TAG", "TCG", "TGG", "TTG", "TAT", "TCT", "TGT", "TTT"
        ]
        self.codon_count = []
        # Matrix used to map codon inputs to index locations in the lists codon_count and codons
        self.matrix = [
            [[0,  1,  2,  3], [4,  5,  6,  7], [8,  9,  10, 11], [12, 13, 14, 15]],     # Matrix A
            [[16, 17, 18, 19], [20, 21, 22, 23], [24, 25, 26, 27], [28, 29, 30, 31]],   # Matrix C
            [[32, 33, 34, 35], [36, 37, 38, 39], [40, 41, 42, 43], [44, 45, 46, 47]],   # Matrix G
            [[48, 49, 50, 51], [52, 53, 54, 55], [56, 57, 58, 59], [60, 61, 62, 63]]    # Matrix T
        ]
        self.setup()

    # ----------------------------------------------------------------------------------------------------------------------
    # Populates the list that tracks the number of codons with zeros
    def setup(self):

        for i in range(64):
            self.codon_count.append(0)

    # ----------------------------------------------------------------------------------------------------------------------
    # Matches the base character with an index in the matrix or returns -1
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

        if (len(codon) != 3):
            print("ERROR: NOT A CODON!")
            return
        
        matrix = self.char_to_index(codon[0])
        row = self.char_to_index(codon[1])
        column = self.char_to_index(codon[2])

        if (matrix == -1 or row == -1 or column == -1):
            print("ERROR: INVALID INPUT!")
            return

        index = self.matrix[matrix][column][row]    # Retreives mapped index in the matrix 
        # print("Index: ", index)
        self.codon_count[index] += 1                # Incriments corresponding codon count
            

    # ----------------------------------------------------------------------------------------------------------------------
    # Returns a list of the codons and their count
    def getCount(self):
        return [self.codons, self.codon_count]

    # ----------------------------------------------------------------------------------------------------------------------
    # Tests if elements are being populated properly
    def debugging(self):
        for i in range(64):
            print(self.codons[i], self.codon_count[i], sep=" : ")


if __name__ == "__main__":

    test = map_codons()
    codons = [ # list of all 64 possible DNA codon combinations
        "AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT",
        "CAA", "CCA", "CGA", "CTA", "CAC", "CCC", "CGC", "CTC", "CAG", "CCG", "CGG", "CTG", "CAT", "CCT", "CGT", "CTT",
        "GAA", "GCA", "GGA", "GTA", "GAC", "GCC", "GGC", "GTC", "GAG", "GCG", "GGG", "GTG", "GAT", "GCT", "GGT", "GTT",
        "TAA", "TCA", "TGA", "TTA", "TAC", "TCC", "TGC", "TTC", "TAG", "TCG", "TGG", "TTG", "TAT", "TCT", "TGT", "TTT"
    ]

    for x in codons:
        test.add_codon(x)

    # test.debugging()
    mycodons = test.getCount()

    for x in range(64):
        print(mycodons[0][x], " : ", mycodons[1][x])

    test.add_codon("TTTT")
    test.add_codon("BBB")