# ==========================================================================================================================
# Lookup Table for codon to amino acids
# ==========================================================================================================================
c_to_a = {
            "AAA": "Lys", "AAC": "Asn", "AAG": "Lys", "AAT": "Asn", "ACA": "Thr", "ACC": "Thr", "ACG": "Thr", "ACT": "Thr", "AGA": "Arg", 
            "AGC": "Ser", "AGG": "Arg", "AGT": "Ser", "ATA": "Ile", "ATC": "Ile", "ATG": "Met", "ATT": "Ile", "CAA": "Gln", "CAC": "His", 
            "CAG": "Gln", "CAT": "His", "CCA": "Pro", "CCC": "Pro", "CCG": "Pro", "CCT": "Pro", "CGA": "Arg", "CGC": "Arg", "CGG": "Arg", 
            "CGT": "Arg", "CTA": "Leu", "CTC": "Leu", "CTG": "Leu", "CTT": "Leu", "GAA": "Glu", "GAC": "Asp", "GAG": "Glu", "GAT": "Asp", 
            "GCA": "Ala", "GCC": "Ala", "GCG": "Ala", "GCT": "Ala", "GGA": "Gly", "GGC": "Gly", "GGG": "Gly", "GGT": "Gly", "GTA": "Val", 
            "GTC": "Val", "GTG": "Val", "GTT": "Val", "TAA": "Stp", "TAC": "Tyr", "TAG": "Stp", "TAT": "Tyr", "TCA": "Ser", "TCC": "Ser", 
            "TCG": "Ser", "TCT": "Ser", "TGA": "Stp", "TGC": "Cys", "TGG": "Trp", "TGT": "Cys", "TTA": "Leu", "TTC": "Phe", "TTG": "Leu", 
            "TTT": "Phe"
        }

# ==========================================================================================================================
# Codon class, tracks the codons within a sequence
# ==========================================================================================================================
class Codon:
    def __init__(self, seq:str, header:str):
        # Codon dictionary, track codon count
        self.codon = {
            "AAA": 0, "AAC": 0, "AAG": 0, "AAT": 0, "ACA": 0, "ACC": 0, "ACG": 0, "ACT": 0, "AGA": 0, "AGC": 0, "AGG": 0, "AGT": 0, 
            "ATA": 0, "ATC": 0, "ATG": 0, "ATT": 0, "CAA": 0, "CAC": 0, "CAG": 0, "CAT": 0, "CCA": 0, "CCC": 0, "CCG": 0, "CCT": 0,
            "CGA": 0, "CGC": 0, "CGG": 0, "CGT": 0, "CTA": 0, "CTC": 0, "CTG": 0, "CTT": 0, "GAA": 0, "GAC": 0, "GAG": 0, "GAT": 0, 
            "GCA": 0, "GCC": 0, "GCG": 0, "GCT": 0, "GGA": 0, "GGC": 0, "GGG": 0, "GGT": 0, "GTA": 0, "GTC": 0, "GTG": 0, "GTT": 0,
            "TAA": 0, "TAC": 0, "TAG": 0, "TAT": 0, "TCA": 0, "TCC": 0, "TCG": 0, "TCT": 0, "TGA": 0, "TGC": 0, "TGG": 0, "TGT": 0, 
            "TTA": 0, "TTC": 0, "TTG": 0, "TTT": 0
        }

        self.seq = seq          # Character string of bases
        self.header = header    # Character string for header info
        self.count_codons()     # Begin counting codons

    # ----------------------------------------------------------------------------------------------------------------------
    # Incriments the count of the codon
    def add_codon(self, codon:str):

        if (len(codon) != 3):               # String can only have 3 bases
            print("ERROR: NOT A CODON!")
            return

        self.codon[codon] += 1              # Incriments the codon count

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the genome sequence string into individual codon sub-strings and adds them to the count
    def count_codons(self):
        end = len(self.seq) - len(self.seq) % 3   # Ignore any trailing characters

        for i in range(0, end, 3):      # Iterate thru each 3 chars in the string
            codon = self.seq[i:(i+3)]
            self.add_codon(codon)       # Add codon to count

    # ----------------------------------------------------------------------------------------------------------------------
    # Returns the totals of all codons as a list [codon, count]
    def get_codon_count(self):
        data = []
        
        for i in self.codon.keys():
            data.append([i, self.codon[i]])
    
        return data

    # ----------------------------------------------------------------------------------------------------------------------
    # Testing if elements are being populated properly
    def debugging(self):
        for i in self.codon:
            print(self.codon[i])

# ==========================================================================================================================
# Amino_Acid class, handles sequence data manipulation
# ==========================================================================================================================
class Amino_Acid:
    def __init__(self, seq:str, header:str):

        # Amino acid dictionary, tracks amino acid count
        self.amino_acid = { 
            "Leu": {"letter": "L", "count": 0}, "Val": {"letter": "V", "count": 0}, "Thr": {"letter": "T", "count": 0}, "Ala": {"letter": "A", "count": 0}, 
            "Ser": {"letter": "S", "count": 0}, "Gly": {"letter": "G", "count": 0}, "Lys": {"letter": "K", "count": 0}, "Asn": {"letter": "N", "count": 0}, 
            "Ile": {"letter": "I", "count": 0}, "Asp": {"letter": "D", "count": 0}, "Phe": {"letter": "F", "count": 0}, "Glu": {"letter": "E", "count": 0}, 
            "Tyr": {"letter": "Y", "count": 0}, "Pro": {"letter": "P", "count": 0}, "Gln": {"letter": "Q", "count": 0}, "Arg": {"letter": "R", "count": 0}, 
            "Cys": {"letter": "C", "count": 0}, "Met": {"letter": "M", "count": 0}, "His": {"letter": "H", "count": 0}, "Trp": {"letter": "W", "count": 0}, 
            "Stp": {"letter": "O", "count": 0}
        }

        self.seq = seq
        self.header = header
        self.count_amino_acids()

    # ----------------------------------------------------------------------------------------------------------------------
    # Incriments the count of the amino acid
    def add_amino_acid(self, amino_acid:str):

        if (len(amino_acid) != 3):                  # String can only have 3 bases
            print("ERROR: NOT A CODON!")
            return

        self.amino_acid[amino_acid]["count"] += 1   # Incriments the amino acid count

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the genome sequence string into individual codon sub-strings, converts them to acids, and adds them to the count
    def count_amino_acids(self):
        end = len(self.seq) - len(self.seq) % 3     # Ignore any trailing characters

        for i in range(0, end, 3):                  # Iterate thru each 3 chars in the string
            codon = self.seq[i:(i+3)]
            self.add_amino_acid(c_to_a[codon])      # Add amino acid to count

    # ----------------------------------------------------------------------------------------------------------------------
    # Converts a Codon sequence to an amino acid sequence and returns it
    def codon_to_amino(self, stops:bool=False):
        amino_str = ""

        start = self.seq.find("ATG")                            # Find index of first start codon
        if(not stops or start < 0):                             # Start at the beginning of the sequence
            start = 0

        #print("I: ", start, " J: ", stop)
        for i in range(start, len(self.seq), 3):                # Convert selected codons to amino acid sequences
            codon = self.seq[i:(i+3)]                           # Get codon
            amino = c_to_a[codon]                               # Convert codon to amino acid

            if(stops and amino == "Stp"):                       # Stop codon reached (TAA, TAG, TGA)
                amino_str += self.amino_acid[amino]["letter"]
                break

            amino_str += self.amino_acid[amino]["letter"]       # Append amino acid letter to string

        return amino_str
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Returns the totals of all amino acids as a list [amino acid, count]
    def get_amino_count(self):
        amino_acid = []

        for key in self.amino_acid.keys():
            amino_acid.append([key, self.amino_acid[key]["count"]])

        return amino_acid
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Testing if elements are being populated properly
    def debugging(self):
        for i in self.amino_acid:
            print(self.amino_acid[i])

# ==========================================================================================================================
# Kmerclass, handles k-mer data manipulation
# ==========================================================================================================================
class Kmer:
    def __init__(self, seq:str, header:str, k:int):
        self.seq = seq          # Input sequence
        self.header = header    # Sequence header information
        self.k = k              # Size of the k-mer (string length)

        self.kmers = {}         # Dictionary that stores the kmer and the count
        self.count_kmers()

    # ----------------------------------------------------------------------------------------------------------------------
    # Adds the k-mer to the count
    def add_kmer(self, kmer:str):

        if(len(kmer) != self.k):    # K-mer must have the same length k
            print("ERROR: K-mer has incorrect length!")
            return

        if kmer in self.kmers:      # K-mer has already been added
            self.kmers[kmer] += 1
        else:                       # K-mer needs to be added
            self.kmers[kmer] = 1

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the k-mer from the sequence and adds it to the count
    def count_kmers(self):
        end = (len(self.seq) - len(self.seq)) % self.k  # Ignore any trailing characters

        for i in range(0, end, 3):                      # Iterate thru each 3 chars in the string
            kmer = self.seq[i:(i+self.k)]
            self.add_kmer(kmer)                         # Add amino acid to count

    # ----------------------------------------------------------------------------------------------------------------------
    # Returns the totals of all k-mers as a list [k-mer, count]
    def get_kmer_count(self):
        kmers = []

        for key in self.kmers.keys():
            kmers.append([key, self.kmers[key]])

        return kmers
    
# ==========================================================================================================================
def main():

    test = Codon()
    codons = [ # list of all 64 possible DNA codon combinations
        "AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT",
        "CAA", "CCA", "CGA", "CTA", "CAC", "CCC", "CGC", "CTC", "CAG", "CCG", "CGG", "CTG", "CAT", "CCT", "CGT", "CTT",
        "GAA", "GCA", "GGA", "GTA", "GAC", "GCC", "GGC", "GTC", "GAG", "GCG", "GGG", "GTG", "GAT", "GCT", "GGT", "GTT",
        "TAA", "TCA", "TGA", "TTA", "TAC", "TCC", "TGC", "TTC", "TAG", "TCG", "TGG", "TTG", "TAT", "TCT", "TGT", "TTT"
    ]

    for x in codons:
        test.add_codon(x)
        #test.add_to_sequence(x)
    
    t = Codon()
    t.add_to_sequence("AATGAAAAAAAAAAAATAAA")
    # print(t.get_amino_count())
    print(t.codon_to_amino())

if __name__ == "__main__":
    main()
    