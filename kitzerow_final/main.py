# Final project

import sys
import numpy as np
import de_bruijn as db

# ==============================================================================================================
# Retreives header/sequence pair data from the fna file
# Retruns seq_data: list of sequences, header_data: List of corresponding header info
def get_data(filename):

    seq_data = []
    header_data = []
    mime = filename.split('.').pop()                        # Get file MIME type

    if(mime == 'fna'):                                      # Input file is a fna (data every two lines)
        print("Reading FNA file: ", filename)

        with open(filename) as f:
            for head, seq in zip(f,f):                      # Get header and sequence info
                head, seq = head.strip(), seq.strip()       # Strip newline characters
                header_data.append(head[1:])                # Add to header list
                seq_data.append(seq)                        # Add to sequece list

    elif(mime == 'fastq'):                                  # File is a fastq (data every four lines)
        print("Reading FASTQ file: ", filename)

        with open(filename) as f:
            for head, seq, p, score in zip(f,f,f,f):        # Get four line at a time (Header, sequence, plus thingy, score)
                head, seq = head.strip(), seq.strip()       # Strip header and sequece data of newline chars
                header_data.append(head[1:])                # Add to header list
                seq_data.append(seq)                        # Add to sequence list

    else:
        print("ERROR: Invalid File Type!")
        print("Only fna or fastq types, entered file type is ", mime)
        return

    return (seq_data, header_data)

# ==============================================================================================================
# Writes data to a text file
def make_txt(data, filename='./output/output.txt'):

    with open(filename, 'w', newline='') as file:
        for x in data:
            file.write(str(x) + '\n')

# ==============================================================================================================
def loop_kmer(data, start, end):
    db_data = db.De_bruijn(data[0], data[1])                                       # Create de Bruijn graph

    for i in range(start, end+1):
        db_data.de_bruijn_graph(k=i)
        db_data.make_docs(True, True, True, str(i))
    

# ==============================================================================================================
def main():
    
    # ----------------------------------------------------------------------------------------------------------
    # User enters prefered file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                                                          # Extracting fastq file name from arguments
        data = get_data(fna_file)                                                       # Processing genome data from fna file
        db_graph = db.De_bruijn(data[0], data[1])                                       # Create de Bruijn graph
        db_graph.matplot_graph(False,True)
    
    # ----------------------------------------------------------------------------------------------------------
    # Use default file
    elif(len(sys.argv) == 1):
        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                                       # Processing genome data from fna file
        k = 10

        db_graph = db.De_bruijn(data[0], data[1], k=k)                                       # Create de Bruijn graph
        db_graph.make_docs(True,True,True, str(k))

    # ----------------------------------------------------------------------------------------------------------
    # Use default file and run kmer loop
    # python .\main start stop
    elif(len(sys.argv) == 3):
        start = sys.argv[1]                                                             # Start at this k-mer
        end = sys.argv[2]                                                               # End at this k-mer

        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                                       # Processing genome data from fna file
        loop_kmer(data, int(start), int(end))                                                           # Loop thru k-mers from start to end


if __name__ == "__main__":
    main()
