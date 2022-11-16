# Samuel Kitzerow, kitze012
# Homework 4, Finding variable genomic regions

import sys
import tree
import phylogeny as py
import numpy as np

# ==============================================================================================================
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

# ==============================================================================================================
# Write distance matrix to a text file
def write_dmatrix(data, filename="./output/genetic-distances.txt"):

    with open(filename, 'w', newline='') as file:
        file.write('\t'.join(data.header) + '\n')               # Write column names from header

        for x in range(len(data.header)):                       # Iterate thru rows in the distance matrix
            s = '\t'.join([str(i) for i in data.dmatrix[x]])    # Join the row data together
            file.write(data.header[x] + '\t' + s + '\n')        # Write the header + data for each corresponding row

# ==============================================================================================================
# Write tree matrix to a text file
def write_data(data, header, filename="./output/edges.txt"):
    tr = tree.Ptree(data, header)   # Initializing tree
    tr.write_edges()                # Writing edges to file
    tr.write_tree()                 # Writing newick tree to file

# ==============================================================================================================
# Main function that calls Pylogeny and write functions
def tree_stuff(header, data):
    data = py.Phylogeny(data, header)            # Initializing neighbor joining
    write_dmatrix(data)
    data.neighbor_joining()
    write_data(data.edges, data.header)

# ==============================================================================================================
def main():
    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from fna file
    # python .\phylogeny.py input_file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                          # Extracting fna file name from arguments
        genes = get_data(fna_file)                      # Processing genome data from fna file
                

    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from example fna file
    # python .\phylogeny.py
    elif(len(sys.argv) == 1):
        fna_file = "./input/Homework4-seqs-with-primers.fna"
        genes = get_data(fna_file)
        

    # ----------------------------------------------------------------------------------------------------------
    else:
        print("ERROR: INVALID ARGUMENMTS!")

# ================================================================================================================================================
if __name__ == "__main__":
    main()