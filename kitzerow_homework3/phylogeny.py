# Samuel Kitzerow, kitze012
# Homework 3, Phylogeny Inference

# Generate the distances between sequences using neighbor joining and build a phylogenic tree

'''
Input: fna data

Output:
1) genetic-distances.txt: A tab-delimited table of pairwise distances between all sequences.
2) edges.txt: A tab delimited. Each row describes an edge in the tree.
3) tree.txt: NEWICK format with all edge distances and with only tips named. For example, (A:0.1,B:0.2,(C:0.3,D:0.4):0.5)

'''
import csv
import sys
import numpy as np

class Phylogeny:
    def __init__(self, seq = [], header = []):
        self.seq = seq
        self.header = header

    def debug(self):
        print(self.seq, self.header)



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
        elif(line[0] == ">"):        # Header information
            header_data.append(line)
        else:                                           # Genetic sequence
            seq_data.append(line)
    
    file.close()
    return (seq_data, header_data)

# =======================================================================================

def main():
    # ----------------------------------------------------------------------------------------------------------
    # Processes codon data from fna file and writes codons with their respective counts to a csv file
    # python .\main.py -c input_file output_file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                              # Extracting fna file name from arguments
        genes = get_data(fna_file)                    # Processing genome data from fna file

    # ----------------------------------------------------------------------------------------------------------
    elif(len(sys.argv) == 2 and sys.argv[1] == "test"):
        fna_file = "hw3.fna"
        genes = get_data(fna_file)

    else:
        print("ERROR: INVALID ARGUMENMTS!")

# ================================================================================================================================================
if __name__ == "__main__":
    main()
    