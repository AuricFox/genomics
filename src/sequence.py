from typing import List
# ==========================================================================================================================
# Lookup Table for codon to amino acids and amino acids to letters
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

a_to_l = { 
            "Leu": "L", "Val": "V", "Thr": "T", "Ala": "A", "Ser": "S", "Gly": "G", "Lys": "K", "Asn": "N", 
            "Ile": "I", "Asp": "D", "Phe": "F", "Glu": "E", "Tyr": "Y", "Pro": "P", "Gln": "Q", "Arg": "R", 
            "Cys": "C", "Met": "M", "His": "H", "Trp": "W", "Stp": "O"
        }

# ==========================================================================================================================
# Codon class, tracks the codons within a sequence
# ==========================================================================================================================
class Codon:
    """
    Parameter(s):
        * header (str): information detailing the genetic sequence
        * seq (List[str]): genetic sequences being evaluated for counting
    """
    def __init__(self, seq: List[str], header:str):
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
    """
    Incriments the count of the codon
    Parameter(s):
        * codon (str): codon being addded to the total count
    """
    def add_codon(self, codon:str):

        if (len(codon) != 3):               # String can only have 3 bases
            print("ERROR: CODON HAS INCORRECT LENGTH!")
            return
        
        if codon not in self.codon:         # Bases can only be A C G or T
            print("ERROR: NOT A CODON!")
            return

        self.codon[codon] += 1              # Incriments the codon count

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the genome sequence string into individual codon sub-strings and adds them to the count
    def count_codons(self):
        for sequence in self.seq:                       # Loop thru all the sequences
            end = len(sequence) - len(sequence) % 3     # Ignore any trailing characters

            for i in range(0, end, 3):                  # Iterate thru each 3 chars in the string
                codon = sequence[i:(i+3)]
                self.add_codon(codon)                   # Add codon to count

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Parses the genome sequence string into individual codon sub-strings and adds them to the count
    Parameter(s):
        * seq (str): genetic sequence being parse into individual codons for counting
    """
    def count_codons_str(self, seq:str):
        end = len(seq) - len(seq) % 3       # Ignore any trailing characters

        for i in range(0, end, 3):          # Iterate thru each 3 chars in the string
            codon = seq[i:(i+3)]
            self.add_codon(codon)           # Add codon to count

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

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Adds the codon totals from another codon object to self's codon totals
    Parameter(s):
        * other (Codon): contains the codon totals of another Codon object
    """
    def __add__(self, other: 'Codon'):
        for key, value in other.codon.items():
            self.codon[key] = self.codon.get(key, 0) + value

# ==========================================================================================================================
# Amino_Acid class, handles sequence data manipulation
# ==========================================================================================================================
class Amino_Acid:
    """
    Parameter(s):
        * seq (List[str]): genetic sequences being evaluated for counting
        * header (str): information detailing the genetic sequence
    """
    def __init__(self, seq: List[str], header:str):

        # Amino acid dictionary, tracks amino acid count
        self.amino_acid = { 
            "Leu": 0, "Val": 0, "Thr": 0, "Ala": 0, "Ser": 0, "Gly": 0, "Lys": 0, "Asn": 0, "Ile": 0, "Asp": 0, "Phe": 0, "Glu": 0, 
            "Tyr": 0, "Pro": 0, "Gln": 0, "Arg": 0, "Cys": 0, "Met": 0, "His": 0, "Trp": 0, "Stp": 0
        }

        self.seq = seq
        self.header = header
        self.count_amino_acids()

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Incriments the count of the amino acid
    Parameter(s):
        * amino_acid (str): amino acid being added to the total count
    """
    def add_amino_acid(self, codon:str):

        if (len(codon) != 3):                       # String can only have 3 bases
            print("ERROR: CODON HAS INCORRECT LENGTH!")
            return
        
        if codon not in c_to_a:                     # Bases can only be A C G or T
            print("ERROR: NOT A CODON")
            return

        amino_acid = c_to_a[codon]
        self.amino_acid[amino_acid] += 1   # Incriments the amino acid count

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the genome sequence string into individual codon sub-strings and adds them to the amino acid count
    def count_amino_acids(self):
        for sequence in self.seq:                       # Iterate thru all sequences
            end = len(sequence) - len(sequence) % 3     # Ignore any trailing characters

            for i in range(0, end, 3):                  # Iterate thru each 3 chars in the string
                codon = sequence[i:(i+3)]
                self.add_amino_acid(codon)              # Add amino acid to count

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Parses the genome sequence string into individual codon sub-strings, converts them to acids, and adds them to the count
    Parameter(s):
        * seq(str): genetic sequence being parse into individual codons and converted to amino acids for counting
    """
    def count_amino_acids_str(self, seq:str):
        end = len(seq) - len(seq) % 3               # Ignore any trailing characters

        for i in range(0, end, 3):                  # Iterate thru each 3 chars in the string
            codon = seq[i:(i+3)]
            self.add_amino_acid(codon)              # Add amino acid to count

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Converts a Codon sequence to an amino acid sequence and return it
    Parameter(s):
        * stops(bool): used to for gene expression, use start and stop codons
    """
    def codon_to_amino(self, stops:bool=False):
        amino_str = ""

        start = self.seq.find("ATG")                    # Find index of first start codon
        if(not stops or start < 0):                     # Start at the beginning of the sequence
            start = 0

        #print("I: ", start, " J: ", stop)
        for i in range(start, len(self.seq), 3):        # Convert selected codons to amino acid sequences
            codon = self.seq[i:(i+3)]                   # Get codon
            amino = c_to_a[codon]                       # Convert codon to amino acid

            if(stops and amino == "Stp"):               # Stop codon reached (TAA, TAG, TGA)
                amino_str += a_to_l[amino]
                break

            amino_str += a_to_l[amino]                  # Append amino acid letter to string

        return amino_str
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Returns the totals of all amino acids as a list [amino acid, count]
    def get_amino_count(self):
        amino_acid = []

        for key in self.amino_acid.keys():
            amino_acid.append([key, self.amino_acid[key]])

        return amino_acid
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Testing if elements are being populated properly
    def debugging(self):
        for i in self.amino_acid:
            print(self.amino_acid[i])

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Adds the amino acid totals from another Amino_Acid object to self's amino acid totals
    Parameter(s):
        * other (Amino_Acid): contains the amino acid totals of another Amino_Acid object
    """
    def __add__(self, other: 'Amino_Acid'):
        for key, value in other.amino_acid.items():
            self.amino_acid[key] = self.amino_acid.get(key, 0) + value

# ==========================================================================================================================
# Kmerclass, handles k-mer data manipulation
# ==========================================================================================================================
class Kmer:
    """
    Parameter(s):
        * seq (List[str]): genetic sequences being evaluated for counting
        * header (str): information detailing the genetic sequence
        * k (int): size of the k-mer, length of the string after parsing
    """
    def __init__(self, seq: List[str], header:str, k:int):
        self.seq = seq
        self.header = header
        self.k = k

        self.kmers = {}         # Dictionary that stores the kmer and the count
        self.count_kmers()

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Adds the k-mer to the count
    Parameter(s):
        * kmer (str): k-mer being added to the total count
    """
    def add_kmer(self, kmer:str):

        if(len(kmer) != self.k):                        # K-mer must have the same length k
            print("ERROR: K-mer has incorrect length!")
            return
        
        bases = ['A', 'C', 'G', 'T']
        if not all(base in bases for base in kmer):    # Perform sanitation check
            print("ERROR: BASES CAN ONLY BE A, C, G, or T!", kmer)
            return

        if kmer in self.kmers:                          # K-mer has already been added
            self.kmers[kmer] += 1
        else:                                           # K-mer needs to be added
            self.kmers[kmer] = 1

    # ----------------------------------------------------------------------------------------------------------------------
    # Parses the k-mer from the sequence and adds it to the count
    def count_kmers(self):
        for sequence in self.seq:                           # Iterate thru all sequences
            end = len(sequence) - len(sequence) % self.k    # Ignore any trailing characters

            for i in range(0, end, self.k):                 # Iterate thru each k chars in the string
                kmer = sequence[i:(i+self.k)]
                self.add_kmer(kmer)                         # Add amino acid to count

    # ----------------------------------------------------------------------------------------------------------------------
    """
    Parses the k-mer from the sequence and adds it to the count
    Parameter(s):
        * seq (str): genetic sequence being parse into individual k-mers for counting
    """
    def count_kmer_str(self, seq:str):
        end = len(seq) - len(seq) % self.k          # Ignore any trailing characters

        for i in range(0, end, self.k):             # Iterate thru each k chars in the string
            kmer = seq[i:(i+self.k)]
            self.add_kmer(kmer)                     # Add amino acid to count

    # ----------------------------------------------------------------------------------------------------------------------
    # Returns the totals of all k-mers as a list [k-mer, count]
    def get_kmer_count(self):
        kmers = []

        for key in self.kmers.keys():
            kmers.append([key, self.kmers[key]])

        return kmers
    
    # ----------------------------------------------------------------------------------------------------------------------
    """
    Adds the k-mer totals from another Kmer object to the self's k-mer totals
    Parameter(s):
        * other(Kmer): contains the k-mer totals of another Kmer object
    """
    def __add__(self, other: 'Kmer'):
        for key, value in other.kmers.items():
            self.kmers[key] = self.kmers.get(key, 0) + value
    
# ==========================================================================================================================
def main():

    test = Kmer(['AAAAAAA'], 'KMER', 4)
    print(test.kmers)

if __name__ == "__main__":
    main()
    