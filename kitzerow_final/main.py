# Final project

'''
Process Frow:
De Bruijn Garph -> Eulerian Cycle/Path -> Alignment
     (1)        ->         (2)         ->    (3)

(1) reads.FASTQ -> Graph.pdf, edges.txt, directed_graph.txt
        Converts reads into kmers and edges. This data is then used to constuct a directed graph file and figure.

(2) directed_graph.txt -> eulerianPath.txt, eulerianCycle.txt (assembled contigs)
        Converts a directed graph file into an Eulerian cycle/path file. This is an assembled contig.

(3) eulerianPath.txt, eulerianCycle.txt -> align.txt, comparison.pdf
        Compares the Eulerian cycle/path file with the assembled spike protein and creates text file with the comparisons and a plot.
'''

import sys
import os
import numpy as np
import de_bruijn as db
import alignment as al
import eulerianCycle as ec
import eulerianPath as ep

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

    elif(mime == 'txt'):                                    # File is a Text
        print("Reading text file ", filename)

        with open(filename) as f:
            for seq in f:                                   # Read each line
                seq = seq.strip()                           # Strip data of new line characters
                seq_data.append(seq)                        # Append data to list

            header_data.append('Assembled data')            # No headers should be in the file so add this one
    else:
        print("ERROR: Invalid File Type!")
        print("Only fna, fastq, or txt types! The entered file type is ", mime)
        return

    return (seq_data, header_data)

# ==============================================================================================================
# Writes data to a text file
def make_txt(data, filename='./output/temp/output.txt'):

    with open(filename, 'w', newline='') as file:
        for x in data:
            file.write(str(x) + '\n')

# ==============================================================================================================
# Loops thru a k-mer range, builds a set of kmers and edges
def loop_kmer(data, kstart, kend, lstart=0, lend=1):
    res = input("\nDo you want alignment files (y)?: ")
    db_data = db.De_bruijn(data[0], data[1])                    # Create de Bruijn graph

    for i in range(kstart, kend+1):                             # Loop thru k-mer range
        db_data.de_bruijn_graph(start=lstart, end=lend, k=i)    # Create de Bruijn graph
        file = db_data.make_docs(True, True, True, str(i))      # Create all files and get directed file

        if(res == 'Y' or res == 'y'):                           # Perform alignment and create corresponding files
            files = eulerian_string(file, str(i))               # Create Eulerian cycle and path
            align_contig(files[0], i)

# ==============================================================================================================
# Clears out old files
def rmove():
    dir = './output/align/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in align directory

    dir = './output/eulerian/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in eulerian directory

    dir = './output/graph/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in graph directory

    dir = './output/temp/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in temp directory

# ==============================================================================================================
# Finds eulerian cycle and path, writes the data to a text file
def eulerian_string(filename, kmer):
    print("Reading: " + filename)

    text = []
    with open(filename, 'r') as f:                              # Open directed graph file
        for line in f:
            text.append(line.strip())                           # Add the edges and strip the newline characters

    cycle = ec.EulerianCycle(text)                              # Find eulerian cycle
    path = ep.EulerianPath(text)                                # Find eulerian cycle
    cycle_file = 'output/eulerian/eulerianCycle_' + kmer + '.txt'  # File that will be written to
    path_file = 'output/eulerian/eulerianPath_' + kmer + '.txt'    # File that will be written to

    with open(cycle_file, 'w') as f:                            # Write eulerian cycle to file

        contig = cycle[0]
        for c in cycle[1:]:                                     # Assembling a contig
            contig += c[len(c)-1]                               # Get last character and add it to contig

        f.write(''.join(contig))

    with open(path_file, 'w') as f:                             # Write eulerian path to file

        contig = path[0]
        for c in path[1:]:                                      # Assembling a contig
            contig += c[len(c)-1]                               # Get last character and add it to contig

        f.write(''.join(path))

    return (cycle_file, path_file)

# ==============================================================================================================
# Aligns assembled contig with reference genome, creates text file of alignment, and plots comparison
def align_contig(filename, kmer):
    ref_seg = get_data('./input/sars_spike_protein_assembled.fna')      # Get data for reference sequence (Spike protien)
    data = get_data(filename)                                           # Get in data from files for alignment

    a = al.alignment(ref_seg[0][0], data[0][0], -2, -1, True)                 # Initialize alignment
    a.alignment_file(kmer)
    a.plot_compare(kmer)

# ==============================================================================================================
def main():

    # ----------------------------------------------------------------------------------------------------------
    # User enters prefered file
    if(len(sys.argv) == 2):
        rmove()
        inpt = sys.argv[1]                                          # Extracting fastq file name from arguments
        file = "./input/sars_spike_protein_reads.fastq"

        if(inpt == '-l'):                                           # User enters one kmer value and read range
            k = int(input("Enter k-mer value: "))
            lstart = int(input("Enter index for first sequence: "))
            lend = int(input("End index for the last sequence: "))

            data = get_data(file)                                   # Processing genome data from fna file
            loop_kmer(data, k, k, lstart, lend)

        elif(inpt == '-k'):                                           # User enters kmer range
            kstart = int(input("Enter starting k-mer value: "))
            kend = int(input("Enter ending k-mer value: "))

            data = get_data(file)                                   # Processing genome data from fna file
            loop_kmer(data, kstart, kend)

        elif(inpt == '-kl'):                                        # User enters kmer range and read range
            kstart = int(input("Enter starting k-mer value: "))
            kend = int(input("Enter ending k-mer value: "))
            lstart = int(input("Enter index for first sequence: "))
            lend = int(input("End index for the last sequence: "))

            data = get_data(file)                                   # Processing genome data from fna file
            loop_kmer(data, kstart, kend, lstart, lend)

        else:
            print("ERROR: Invalid Input!")
            print("Use:\n-l\tsegmenting reads list\n-k\tk-mer range\n-kl\tk-mer range and segmenting reads list\n-a\taligning contig")
            return
    # ----------------------------------------------------------------------------------------------------------
    # Use default file
    elif(len(sys.argv) == 1):
        rmove()

        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                       # Processing genome data from fna file
        k = 14
        lstart = 0
        lend = 4
        
        loop_kmer(data, k, k, lstart, lend)

    # ----------------------------------------------------------------------------------------------------------
    # Use default file and run kmer loop
    # python .\main start stop
    '''
    elif(len(sys.argv) == 3):
        rmove()

        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                       # Processing genome data from fna file

    '''
if __name__ == "__main__":
    main()
