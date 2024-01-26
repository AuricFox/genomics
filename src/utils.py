import os, re, mimetypes, json, csv, zipfile, logging
from typing import List

PATH = os.path.dirname(os.path.abspath(__file__))

logging.basicConfig(
    filename=os.path.join(PATH, './output/app.log'),
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s]: %(message)s'
)

LOGGER = logging.getLogger(__name__)

# ========================================================================================================================================
# Functions used for processing files
# ========================================================================================================================================
def create_file(file):
    '''
    Takes in a file object, sanitizes, validates, and saves it to the temp directory

    Parameter(s):
        file: the user input file being saved

    Output(s):
        file_path (str): the path to the saved file
    '''

    # Replace special characters with underscores
    sanitized_name = re.sub(r'[\\/*?:"<>|]', '_', file.filename)
    # Remove leading and trailing whitespace
    sanitized_name = sanitized_name.strip()
    # Getting the name of the file without the extension
    sanitized_name = sanitized_name.split('.')[0]

    allowed_mime_types = ['image/jpeg', 'image/png', 'application/pdf', 'text/plain']
    allowed_extensions = ['.jpg', '.jpeg', '.png', '.pdf', '.fna', '.fastq', '.txt']

    # Get the file's MIME type and extension
    file_mime_type, _ = mimetypes.guess_type(file.filename)
    file_extension = os.path.splitext(file.filename)[1].lower()

    # Check if the file's MIME type or extension is allowed
    if file_mime_type is not None and file_mime_type not in allowed_mime_types:
        logging.error(f'{file.filename} MIME type is not supported! MIME type: {file_mime_type}')
        return None

    if file_extension not in allowed_extensions:
        logging.error(f'{file.filename} extension is not supported! Extension: {file_extension}')
        return None
    
    path = os.path.join(PATH, "temp")                           # Path where file will be saved
    os.makedirs(path, exist_ok=True)                            # Create path if it doesn't exist
    
    original_file_path = os.path.join(path, f'{sanitized_name}{file_extension}')
    new_file_path = original_file_path
    
    counter = 1
    # loop thru the files to ensure the saved file does not have the same name as another
    while os.path.exists(new_file_path):
        new_file_path = os.path.join(path, f'{sanitized_name}_{counter}{file_extension}')
        counter += 1
    
    logging.info(f"Creating file: {new_file_path}")
    file.save(new_file_path)

    return new_file_path

# ----------------------------------------------------------------------------------------------------------------------------
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
        LOGGER.info(f"Reading FNA file: {filename}")

        with open(filename) as f:
            for head, seq in zip(f,f):                      # Get header and sequence info
                head, seq = head.strip(), seq.strip()       # Strip newline characters
                headers.append(head[1:])                    # Add to header list
                sequences.append(seq)                       # Add to sequece list

    # Handle fastq files (data every four lines)
    elif(mime == 'fastq'):
        LOGGER.info(f"Reading FASTQ file: {filename}")

        with open(filename) as f:
            for head, seq, p, score in zip(f,f,f,f):        # Get four line at a time (Header, sequence, plus thingy, score)
                head, seq = head.strip(), seq.strip()       # Strip header and sequece data of newline chars
                headers.append(head[1:])                    # Add to header list
                sequences.append(seq)                       # Add to sequence list

    # Handle text files
    elif(mime == 'txt'):
        LOGGER.info(f"Reading text file {filename}")

        with open(filename) as f:
            for seq in f:                                   # Read each line
                seq = seq.strip()                           # Strip data of new line characters
                sequences.append(seq)                       # Append data to list

            headers.append('Assembled data')            # No headers should be in the file so add this one
    else:
        LOGGER.error(f"ERROR: Invalid File Type! Only fna, fastq, or txt types! The entered file type is {mime}")
        return

    return (sequences, headers)

# ----------------------------------------------------------------------------------------------------------------------------
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
        LOGGER.info(f"Merging files: {file1}, {file2}\nMerge Destination:{file3}")

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
        
    except FileNotFoundError:
        LOGGER.error("File not found. Please check the input filenames.")
    except Exception as e:
        LOGGER.error(f"An error occurred: {str(e)}")

# ----------------------------------------------------------------------------------------------------------------------------
def create_zip(files:List[str], zipname:str='./temp/output.zip'):
    '''
    Takes in a list of files and zips them up into one zipfile

    Parameter(s):
        files (List[str]): list of file names being zipped
        zipname (str, default=./temp/output.zip): path of the zip file being saved (path from src)

    Output(s):
        file_path (str): the path to the saved zip file or None if no files to process
    '''
    try:
        file_path = os.path.join(PATH, zipname)             # Creating saved file path
    
        if len(files) == 0:                                 # No files to process
            return None
    
        with zipfile.ZipFile(file_path, "w") as zipf:
            for file in files:
                zipf.write(file, os.path.basename(file))

        LOGGER.info(f"Zipped Files: {files}\nZipped Destination: {file_path}")
        return file_path
    
    except Exception as e:
        LOGGER.error(f"Failed to zip {files}: {str(e)}")

# ========================================================================================================================================
# Functions used for removing files
# ========================================================================================================================================
def remove_file(filename:str):
    '''
    Takes in a file object and removes the file from the temp directory

    Parameter(s):
        file (str): the input file being removed

    Output(s): None
    '''

    try:
        path = os.path.join(os.path.dirname(__file__), "src/temp")  # Path where file is saved
        file_path = os.path.join(path, filename)                    # Creating saved file path
        os.remove(file_path)                                        # File is no longer needed
        LOGGER.info(f"Successfully removed {file_path}")

    except OSError as e:
        LOGGER.error(f'Error while removing {filename}: {str(e)}')

# ----------------------------------------------------------------------------------------------------------------------------
def remove_files(files:List[str]):
    '''
    Remove files from memory

    Parameter(s):
        files (List[str]): a list of files being removed

    Output(s): None
    '''

    for file in files:

        try:
            file_path = os.path.join(PATH, file)                # Creating saved file path
            os.remove(file_path)                                # File is no longer needed
            LOGGER.info(f"Successfully removed {file_path}")

        except OSError as e:
            LOGGER.error(f'Error while removing {file}: {str(e)}')
# ----------------------------------------------------------------------------------------------------------------------------
def clean_temp():
    '''
    Removes all files from the temp directory.

    Parameter(s): None

    Output(s):
        None, removes all files from the temp folder.
    '''
    try:
        folder = os.path.join(PATH, './temp/')

        for file in os.listdir(folder):
            os.remove(os.path.join(PATH, file))
    
    except PermissionError as e:
        LOGGER.error(f"Permission error when removing files: {str(e)}")
    except OSError as e:
        LOGGER.error(f"An error occurred while removing files from temp folder: {str(e)}")
    except Exception as e:
        LOGGER.error(f"An error occurred while removing files from temp folder: {str(e)}")

# ========================================================================================================================================
# Function(s) used for creating txt, csv, and Json file
# ========================================================================================================================================
def make_txt(data:dict, header:List[str]=[], filename:str='./temp/output.txt'):
    '''
    Takes in a dictionary and writes the data to a txt file
    
    Parameter(s):
        data (dict): dictionary containing the data being written to a file
        header (List[str], default=[]): list of header info for the corresponding data
        filename (str, default=./temp/output.txt): file path where the data will be saved
        
    Output(s):
        A file with the user input filename containing the dictionary data
    '''
    try:
        file_path = os.path.join(PATH, filename)             # Creating saved file path
    
        with open(file_path, 'w', newline='') as file:
            if header != []:
                file.write('\t'.join(header) + '\n')
    
            for key,value in data.items():
                file.write(f'{key}\t{value}\n')

        LOGGER.info(f"Successfully saved data to {file_path}")
        return file_path
    
    except FileNotFoundError as e:
        LOGGER.error(f"{filename} not found when saving data!")
    except PermissionError as e:
        LOGGER.error(f"Permission error when saving data to {filename}!")
    except Exception as e:
        LOGGER.error(f"Failed to save data to {file_path}: {str(e)}")

# ----------------------------------------------------------------------------------------------------------------------------
def make_csv(data:dict, header:List[str]=[], filename:str='./temp/output.csv'):
    '''
    Takes in a dictionary and writes the data to a csv file
    
    Parameter(s):
        data (dict): dictionary containing the data being written to a file
        header (List[str], default=[]): list of header info for the corresponding data
        filename (str, default=./temp/output.csv): file path where the data will be saved
        
    Output(s):
        A file with the user input filename containing the dictionary data
    '''

    try:
        file_path = os.path.join(PATH, filename)

        with open(file_path, 'w', newline='') as file:
            cfile = csv.writer(file)

            if header != []:
                cfile.writerow(header)

            for key,value in data.items():
                cfile.writerow([key, value])

        LOGGER.info(f"Successfully saved data to {file_path}")
        return file_path
    
    except FileNotFoundError as e:
        LOGGER.error(f"{filename} not found when saving data!")
    except PermissionError as e:
        LOGGER.error(f"Permission error when saving data to {filename}!")
    except Exception as e:
        LOGGER.error(f"Failed to save data to {file_path}: {str(e)}")

# ----------------------------------------------------------------------------------------------------------------------------
def make_json(data:dict, filename:str='./temp/output.json'):
    '''
    Takes in a dictionary and writes the data to a json file
    
    Parameter(s):
        data (dict): dictionary containing the data being written to a file
        filename (str, default=./temp/output.json): file path where the data will be saved
        
    Output(s):
        A file with the user input filename containing the dictionary data
    '''
    try:
        file_path = os.path.join(PATH, filename)

        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)

        LOGGER.info(f"Successfully saved data to {file_path}")
        return file_path
    
    except FileNotFoundError as e:
        LOGGER.error(f"{filename} not found when saving data!")
    except PermissionError as e:
        LOGGER.error(f"Permission error when saving data to {filename}!")
    except Exception as e:
        LOGGER.error(f"Failed to save data to {file_path}: {str(e)}")


# ========================================================================================================================================
# Function(s) used for testing
# ========================================================================================================================================
def runtime_csv(data, header:List[str]=[], filename:str='./temp/runtime.csv'):
    '''
    Create csv file with compiled runtimes
    
    Parameter(s):
        data: runtime data being written to the csv file
        header (List[str], default=[]): list of header info on the data
        filename (str, default=./temp/runtime.csv): file path where the data will be saved
    
    Output(s):
        A path to the file with the runtime data
    '''

    try:
        file_path = os.path.join(PATH, filename)            # Creating saved file path

        with open(file_path, 'w', newline='') as file:
            cfile = csv.writer(file)

            # Add header if it is not empty
            if header != []:
                cfile.writerow(header)

            cfile.writerows(data)                           # Wrights codon or amino acid data to csv

        LOGGER.info(f"Successfully saved runtime data to {file_path}")
        return file_path
    
    except FileNotFoundError as e:
        LOGGER.error(f"{filename} not found when saving data!")
    except PermissionError as e:
        LOGGER.error(f"Permission error when saving data to {filename}!")
    except Exception as e:
        LOGGER.error(f"Failed to save runtime data to {file_path}: {str(e)}")

# ========================================================================================================================================
# Error Handling
# ========================================================================================================================================
class InvalidFile(Exception):
    pass

class FileNotFound(Exception):
    pass

class InvalidInput(Exception):
    pass

# ========================================================================================================================================
if __name__ == "__main__":

    '''
    merge_files('SARS-CoV-2_separate_genes.csv', 'SARS-CoV-2_whole_genome.csv', 'combined_codons.csv')
    merge_files('separate_amino_acids.csv', 'whole_amino_acids.csv', 'combined_amino_acids.csv')

    if(len(sys.argv) == 5 and sys.argv[1] == "merge"):    # Merging two csv files into one
        merge_files(sys.argv[2], sys.argv[3], sys.argv[4])
    '''
    print(PATH)
    
