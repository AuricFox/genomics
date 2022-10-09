# Samuel Kitzerow, kitze012
# Sequence class, handles sequence data manipulation

# Codon types: A C G T
# Mapping codons to index values to an array/list used for counting
# 4 designator types with 3 total designators in a string (4^3 = 64 possible codon combinations)

class sequence:
    def __init__(self, seq = "", header = ""):
        # Codon dictionary, track codon count
        self.codon = {
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
        # Amino acid dictionary, tracks amino acid count
        self.amino_acid = { 
            "Leu": {"letter": "L", "count": 0}, "Val": {"letter": "V", "count": 0}, "Thr": {"letter": "T", "count": 0}, "Ala": {"letter": "A", "count": 0}, 
            "Ser": {"letter": "S", "count": 0}, "Gly": {"letter": "G", "count": 0}, "Lys": {"letter": "K", "count": 0}, "Asn": {"letter": "N", "count": 0}, 
            "Ile": {"letter": "I", "count": 0}, "Asp": {"letter": "D", "count": 0}, "Phe": {"letter": "F", "count": 0}, "Glu": {"letter": "E", "count": 0}, 
            "Tyr": {"letter": "Y", "count": 0}, "Pro": {"letter": "P", "count": 0}, "Gln": {"letter": "Q", "count": 0}, "Arg": {"letter": "R", "count": 0}, 
            "Cys": {"letter": "C", "count": 0}, "Met": {"letter": "M", "count": 0}, "His": {"letter": "H", "count": 0}, "Trp": {"letter": "W", "count": 0}, 
            "Stp": {"letter": "O", "count": 0}
        }

        self.sequence = seq         # Character string of bases
        self.header = header        # Character string for header info

    # ----------------------------------------------------------------------------------------------------------------------
    # Converts a Codon sequence to an amino acid sequence and returns it
    def codon_to_amino(self):
        amino_str = ""

        start = self.sequence.find("ATG")                                                               # Find index of first start codon
        if(start < 0):                                                                                  # No start codon found, start at beginning
            start = 0

        indices = [self.sequence.find("TAA"), self.sequence.find("TAG"), self.sequence.find("TGA")]     # Find indices of stop codons
        indices = [n for n in indices if n >= 0 and n > start+3]                                        # Get non-negative indicies

        if(indices == []):                                                                              # No stop codon found, use sequence length
            stop = len(self.sequence)

        stop = min(indices) + 3                                                                         # Find index of first end codon

        print(self.sequence)
        print("I: ", start, " J: ", stop)
        for i in range(start, stop, 3):                                                                 # Convert selected codons to amino acid sequences
            codon = self.sequence[i:(i+3)]
            amino = self.codon[codon]["amino_acid"]
            amino_str += self.amino_acid[amino]["letter"]

        return amino_str

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the codon string into individual characters, maps them to the matching codon, and incriments the count of the codon
    def add_codon(self, codon):

        if (len(codon) != 3):                                               # String can only have 3 bases
            print("ERROR: NOT A CODON!")
            return

        # print("Index: ", index)
        self.codon[codon]["count"] += 1                                     # Incriments the codon count
        self.amino_acid[self.codon[codon]["amino_acid"]]["count"] += 1      # Add codon to amino acid

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
    # Returns the totals of all codons as a list [codon, count]
    def get_codon_count(self):
        data = []
        
        for i in self.codon.keys():
            data.append([i, self.codon[i]["count"]])        # [codon, count]
    
        return data

    # ----------------------------------------------------------------------------------------------------------------------
    # Returns the totals of all amino acids as a list [codon, count]
    def get_amino_count(self):
        amino_acid = []

        for key in self.amino_acid.keys():
            amino_acid.append([key, self.amino_acid[key]["count"]])

        return amino_acid

    # ----------------------------------------------------------------------------------------------------------------------
    # Testing if elements are being populated properly
    def debugging(self):
        for i in self.codon:
            print(self.codon[i])

        for i in self.amino_acid:
            print(self.amino_acid[i])

# ==========================================================================================================================
def main():

    test = sequence()
    codons = [ # list of all 64 possible DNA codon combinations
        "AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT",
        "CAA", "CCA", "CGA", "CTA", "CAC", "CCC", "CGC", "CTC", "CAG", "CCG", "CGG", "CTG", "CAT", "CCT", "CGT", "CTT",
        "GAA", "GCA", "GGA", "GTA", "GAC", "GCC", "GGC", "GTC", "GAG", "GCG", "GGG", "GTG", "GAT", "GCT", "GGT", "GTT",
        "TAA", "TCA", "TGA", "TTA", "TAC", "TCC", "TGC", "TTC", "TAG", "TCG", "TGG", "TTG", "TAT", "TCT", "TGT", "TTT"
    ]

    for x in codons:
        test.add_codon(x)
        #test.add_to_sequence(x)
    
    t = sequence()
    t.add_to_sequence("AATGAAAAAAAAAAAATAAA")
    # print(t.get_amino_count())
    print(t.codon_to_amino())

if __name__ == "__main__":
    main()
    