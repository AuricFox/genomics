# Samuel Kitzerow, kitze012
# Homework 3, Phylogeny Inference

'''
Generate the distances between sequences using neighbor joining and build a phylogenic tree.
Gap in both sequences are treated as a similarity.

Input: fna data

Output:
1) genetic-distances.txt: A tab-delimited table of pairwise distances between all sequences.
2) edges.txt: A tab delimited. Each row describes an edge in the tree.
3) tree.txt: NEWICK format with all edge distances and with only tips named. For example, (A:0.1,B:0.2,(C:0.3,D:0.4):0.5)

'''
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize, precision=10)

class Phylogeny:
    def __init__(self, seq = [], header = []):
        self.seq = seq                  # List of taxa sequences
        self.header = header            # Correspoding markers for sequences
        self.dmatrix = None             # Distance matrix
        self.edges = None               # Edges of phylogenic tree
        self.tree = None                # Nodes of phylogenic tree

        self.build_Dmatrix()            # Build distance matrix

    # ----------------------------------------------------------------------------------------------------------
    # Calculates % dissimilarity between two sequences
    def distances(self, seq1, seq2):
        match = 0
        mismatch = 0

        if(len(seq1) == 0 or len(seq2) == 0):           # Sequences cannot be null
            print("ERROR: NULL STRING!")
            return
        if(len(seq1) != len(seq2)):
            print("ERROR: lengths are not equal!")      # Lengths of comparing sequences are not equal
        
        total = len(seq1)                               # Total number of bases
        for (i,j) in zip(seq1,seq2):
            if(i == j):
                match += 1
            else:
                mismatch += 1

        return (mismatch/total)                         # Return % dissimilarity


    # ----------------------------------------------------------------------------------------------------------
    # Builds distance matrix of the input sequences
    def build_Dmatrix(self):
        size = len(self.seq)                                            # Size of n x n distance matrix
        self.dmatrix = np.array([[0.0]*size for i in range(size)])      # Initialize distance matrix

        for i in range(size):                                           # Sequence 1
            for j in range(size):                                       # Sequence 2 (Compare)
                dis = self.distances(self.seq[i],self.seq[j])           # Getting % dissimilarity between sequences
                #print(i , ", ", j, ", ", dis)
                self.dmatrix[i][j] = dis                                # Setting % dissimilarity in distance matrix

    # ----------------------------------------------------------------------------------------------------------
    # Prints out values for debugging
    def debug(self):
        #print(self.seq, self.header)
        print(self.dmatrix)


# =======================================================================================
# Retreives data from the fna file and returns a tuple containing a list of sequences and headers
def get_data(filename):
    file = open(filename, 'r')
    seq_data = []
    header_data = []
    
    for line in file:                                   # Read each line in file
        line = line.strip()                             # Strip newline characters
        
        if(line == ""):                                 # Empty line
            # print("Not text")
            continue
        elif(line[0] == ">"):                           # Header information
            header_data.append(line[1:])
        else:                                           # Genetic sequence
            seq_data.append(line)
    
    file.close()
    return (seq_data, header_data)

# =======================================================================================
# Write distance matrix to text file
def write_dmatrix(data, filename="./output/genetic-distances.txt"):

    with open(filename, 'w', newline='') as file:
        file.write('\t'.join(data.header) + '\n')               # Write column names from header

        for x in range(len(data.header)):                       # Iterate thru rows in the distance matrix
            s = '\t'.join([str(i) for i in data.dmatrix[x]])    # Join the row data together
            file.write(data.header[x] + '\t' + s + '\n')        # Write the header + data for each corresponding row

# =======================================================================================
def main():
    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from fna file
    # python .\phylogeny.py input_file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                          # Extracting fna file name from arguments
        genes = get_data(fna_file)                      # Processing genome data from fna file
        data = Phylogeny(genes[0], genes[1])            # Initializing neighbor joining

    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from example fna file
    # python .\phylogeny.py
    elif(len(sys.argv) == 1):
        fna_file = "./example/example.fna"
        genes = get_data(fna_file)
        data = Phylogeny(genes[0], genes[1])

        #data.debug()
        write_dmatrix(data)

    # ----------------------------------------------------------------------------------------------------------
    else:
        print("ERROR: INVALID ARGUMENMTS!")

# ================================================================================================================================================
if __name__ == "__main__":
    main()
    