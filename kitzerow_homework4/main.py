# Samuel Kitzerow, kitze012
# Homework 4, Finding variable genomic regions

import sys
#import tree
import phylogeny as py
import variance as va
import numpy as np

# ==============================================================================================================
# Retreives header/sequence pair data from the fna file
# Retruns seq_data: list of sequences, header_data: List of corresponding header info
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
# Main function that calls Pylogeny and write functions
def tree_stuff(header, data):
    data = py.Phylogeny(data, header)            # Initializing neighbor joining
    write_dmatrix(data)
    data.neighbor_joining()
    #tree.write_data(data.edges, data.header)

# ==============================================================================================================
def main():
    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from fna file
    # python .\phylogeny.py input_file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                          # Extracting fna file name from arguments
        genes = get_data(fna_file)                      # Processing genome data from fna file

        s = va.Variance(genes[0], genes[1])
        s.debug()
                

    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from example fna file
    # python .\phylogeny.py
    elif(len(sys.argv) == 1):
        fna_file = "./input/Homework4-seqs-with-primers.fna"
        genes = get_data(fna_file)

        s = va.Variance(genes[0], genes[1])
        s.debug()

    # ----------------------------------------------------------------------------------------------------------
    else:
        print("ERROR: INVALID ARGUMENMTS!")

# ================================================================================================================================================
if __name__ == "__main__":
    main()