import sys
import csv
import os
import zipfile
import json
from typing import List

PATH = os.path.dirname(os.path.abspath(__file__))

# ========================================================================================================================================
def get_data(filename:str):
    '''
    Retreives header/sequence pair data from fna, fastq, and txt files
    
    Parameter(s):
        filename (str): file that the header and sequence data is being read from
    
    Output(s):
        sequences, headers (List[str],List(str)): a tuple containing lists of sequences and header info
    '''

    sequences = []
    headers = []
    mime = filename.split('.').pop()                        # Get file MIME type

    # Handle fna files (data every two lines)
    if(mime == 'fna'):
        print("Reading FNA file: ", filename)

        with open(filename) as f:
            for head, seq in zip(f,f):                      # Get header and sequence info
                head, seq = head.strip(), seq.strip()       # Strip newline characters
                headers.append(head[1:])                    # Add to header list
                sequences.append(seq)                       # Add to sequece list

    # Handle fastq files (data every four lines)
    elif(mime == 'fastq'):
        print("Reading FASTQ file: ", filename)

        with open(filename) as f:
            for head, seq, p, score in zip(f,f,f,f):        # Get four line at a time (Header, sequence, plus thingy, score)
                head, seq = head.strip(), seq.strip()       # Strip header and sequece data of newline chars
                headers.append(head[1:])                    # Add to header list
                sequences.append(seq)                       # Add to sequence list

    # Handle text files
    elif(mime == 'txt'):
        print("Reading text file ", filename)

        with open(filename) as f:
            for seq in f:                                   # Read each line
                seq = seq.strip()                           # Strip data of new line characters
                sequences.append(seq)                       # Append data to list

            headers.append('Assembled data')            # No headers should be in the file so add this one
    else:
        print("ERROR: Invalid File Type!")
        print(f"Only fna, fastq, or txt types! The entered file type is {mime}")
        return

    return (sequences, headers)

# ========================================================================================================================================
def merge_files(file1:str, file2:str, file3:str):
    '''
    Merges two files together and writes the results to a third file
    
    Parameter(s):
        file1 (str): path of the first file to be mergered
        file2 (str): path of the second file to be merged
        file3 (str): path of the file that the results are written too
    
    Output(s):
        file3 (str): a file containing the data from file1 and file2
    '''

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
def create_zip(files:List[str], zipname:str='./temp/output.zip'):
    '''
    Takes in a list of files and zips them up into one zipfile

    Parameter(s):
        files (List[str]): list of file names being zipped
        zipname (str, optional): path of the zip file being saved (path from src)

    Output(s):
        file_path (str): the path to the saved zip file or None if no files to process
    '''

    if len(files) == 0:                                 # No files to process
        return None
    
    file_path = os.path.join(PATH, zipname)             # Creating saved file path

    with zipfile.ZipFile(file_path, "w") as zipf:
        for file in files:
            zipf.write(file, os.path.basename(file))

    return file_path

# ========================================================================================================================================
def remove_files(files:List[str]):
    '''
    Remove files from memory

    Parameter(s):
        files (List[str]): a list of files being removed

    Output(s): None
    '''

    for file in files:
        print(f"Removing file {file}")

        file_path = os.path.join(PATH, file)                    # Creating saved file path

        try:
            os.remove(file_path)                                # File is no longer needed
        except OSError as e:
            print(f'Error while removing file {file}: {e}')

# ========================================================================================================================================
def make_txt(data:dict, header:List[str]=[], filename:str='./temp/output.txt'):
    '''
    Takes in a dictionary and writes the data to a txt file
    
    Parameter(s):
        data (dict): dictionary containing the data being written to a file
        header (List[str]): list of header info for the corresponding data
        filename (str): file path where the data will be saved
        
    Output(s):
        A file with the user input filename containing the dictionary data
    '''

    file_path = os.path.join(PATH, filename)             # Creating saved file path

    with open(file_path, 'w', newline='') as file:
        if header != []:
            file.write('\t'.join(header) + '\n')

        for key,value in data.items():
            file.write(f'{key}\t{value}\n')

    return file_path

# ========================================================================================================================================
def make_csv(data:dict, header:List[str]=[], filename:str='./temp/output.txt'):
    '''
    Takes in a dictionary and writes the data to a csv file
    
    Parameter(s):
        data (dict): dictionary containing the data being written to a file
        header (List[str]): list of header info for the corresponding data
        filename (str): file path where the data will be saved
        
    Output(s):
        A file with the user input filename containing the dictionary data
    '''

    file_path = os.path.join(PATH, filename)

    with open(file_path, 'w', newline='') as file:
        cfile = csv.writer(file)

        if header != []:
            cfile.writerow(header)

        for key,value in data.items():
            cfile.writerow([key, value])

    return file_path

# ========================================================================================================================================
def make_json(data:dict, filename:str='./temp/output.json'):
    '''
    Takes in a dictionary and writes the data to a json file
    
    Parameter(s):
        data (dict): dictionary containing the data being written to a file
        filename (str): file path where the data will be saved
        
    Output(s):
        A file with the user input filename containing the dictionary data
    '''

    file_path = os.path.join(PATH, filename)

    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

    return file_path

# ========================================================================================================================================
def runtime_csv(data, header:List[str]=[] ,filename:str='./temp/runtime.csv'):
    '''
    Create csv file with compiled runtimes
    
    Parameter(s):
        data: runtime data being written to the csv file
        header (List[str]): list of header info on the data
        filename (str, optional): file path where the data will be saved
    
    Output(s):
        A file with the filename containing the runtime data
    '''
    file_path = os.path.join(PATH, filename)            # Creating saved file path

    with open(file_path, 'w', newline='') as file:
        cfile = csv.writer(file)

        # Add header if it is not empty
        if header != []:
            cfile.writerow(header)

        cfile.writerows(data)                           # Wrights codon or amino acid data to csv

# ========================================================================================================================================
if __name__ == "__main__":

    '''
    merge_files('SARS-CoV-2_separate_genes.csv', 'SARS-CoV-2_whole_genome.csv', 'combined_codons.csv')
    merge_files('separate_amino_acids.csv', 'whole_amino_acids.csv', 'combined_amino_acids.csv')

    if(len(sys.argv) == 5 and sys.argv[1] == "merge"):    # Merging two csv files into one
        merge_files(sys.argv[2], sys.argv[3], sys.argv[4])
    '''
    print(PATH)
    
