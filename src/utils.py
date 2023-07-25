import sys
import csv
import os

# ========================================================================================================================================
# Retreives header/sequence pair data from the fna file
# Inputs:
#   * filename(str): file that the header and sequence data is being read from
# Retruns:
#   * seq_data, header_data(List[str],List(str)): a tuple containing lists of sequences and header info
def get_data(filename:str):

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

# ========================================================================================================================================
# Merges two files together and writes the results to a third file
# Inputs:
#   * file1(str): path of the first file to be mergered
#   * file2(str): path of the second file to be merged
#   * file3(str): path of the file that the results are written too
def merge_files(file1:str, file2:str, file3:str):
    try:
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            reader1 = csv.reader(f1)
            reader2 = csv.reader(f2)
            data = []

            for row1, row2 in zip(reader1, reader2):
                # Assuming file1 and file2 have a common identifier in the first column
                common_identifier = row1[0]

                # Merge data based on the common identifier
                merged_row = row1 + [row2[1]]

                data.append(merged_row)

            with open(file3, 'w', newline='') as output_file:
                writer = csv.writer(output_file)
                writer.writerows(data)
        print("Files merged successfully.")
    except FileNotFoundError:
        print("File not found. Please check the input filenames.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

# ========================================================================================================================================
# Clears out old files
def rmove():
    print("Removing files\n")
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

# ========================================================================================================================================
# Writes data to a text file
def make_txt(data, filename:str='./output/temp/output.txt'):

    with open(filename, 'w', newline='') as file:
        for x in data:
            file.write(str(x) + '\n')

# ========================================================================================================================================
# Writes data to a text file
def make_csv(data, filename:str='./output/runtime.csv'):

    with open(filename, 'w', newline='') as file:
        cfile = csv.writer(file)
        cfile.writerow(['K-mer', '# of Reads', 'Run-Time(s)'])  # Wrights codon header to csv
        cfile.writerows(data)                                   # Wrights codon or amino acid data to csv

# ========================================================================================================================================
if __name__ == "__main__":

    # merge_files('SARS-CoV-2_separate_genes.csv', 'SARS-CoV-2_whole_genome.csv', 'combined_codons.csv')
    # merge_files('separate_amino_acids.csv', 'whole_amino_acids.csv', 'combined_amino_acids.csv')

    if(len(sys.argv) == 5 and sys.argv[1] == "merge"):    # Merging two csv files into one
        merge_files(sys.argv[2], sys.argv[3], sys.argv[4])

    
