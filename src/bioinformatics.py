# GUI terminal
# Runs main program that calls other scripts

import sys
import os
import sequence as sq
import alignment as al
import de_bruijn as db
import variance as vr
import utils

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp")

# ==============================================================================================================
def get_sequence_data(
    file:str, 
    codon:bool=True, 
    amino:bool=False, 
    kmer:bool=False, 
    k:int=3, 
    file_type:str='txt'):
    '''
    Extracts the header and sequence data from the input file and counts the totals in the sequence. The data is
    then saved to individual files depending on the file type and zipped into one zip file. The path to the zip
    file is then returned. 
    
    Parameter(s):
        file (str): path to the file containing the sequence data
        codon (bool, default=True): count the number of codons in a sequence if true else don't count them
        amino (bool, default=False): count the number of amino acids in a sequence if true else don't count them
        kmer (bool, default=False): count the number of k-mers in a sequence if true else don't count them
        k (int, default=3): size of the k-mer, number of bases/characters
        file_type (str, defualt=txt): type of file(s) to be returned in the zip file
    
    Output(s):
        Creates a zip file containing a series of files depending on file_type which contain the counts of codons, 
        amino acids, and/or k-mers if no errors occur. Returns a dictionary response containing the path of the zip 
        file or any errors.
    '''
    response = {}
    files = []

    try:
        sequences = utils.get_data(os.path.join(PATH, file))
        data = {}

        # Counting codons
        if codon:
            seq = sq.Codon(sequences=sequences[0], header="Codons")
            data['codons'] = seq.codon
        # Counting amino acids
        if amino:
            seq = sq.Amino_Acid(sequences=sequences[0], header="Amino Acids")
            data['amino acids'] = seq.amino_acid
        # Counting k-mers
        if kmer:
            seq = sq.Kmer(sequences=sequences[0], header="K-mers", k=k)
            data['kmers'] = seq.kmers

        # Create a text file for each sequence type
        if file_type == 'txt':
            for key,value in data.items():
                filename = f'./temp/{key}_count.txt'
                files.append(utils.make_txt(data=value, header=[key, 'Count'], filename=filename))

        # Create a csv file for each sequence type
        elif file_type == 'csv':
            for key,value in data.items():
                filename = f'./temp/{key}_count.csv'
                files.append(utils.make_csv(data=value, header=[key, 'Count'], filename=filename))

        # Create a json file for the sequence data
        elif file_type == 'json':
            filename = f'./temp/sequenece_count.json'
            files.append(utils.make_json(data=data, filename=filename))

        # Zipping the collection of files
        response['zip_file'] = utils.create_zip(files=files, zipname='./temp/sequence.zip')

    except Exception as e:
        response['error'] = str(e)
    finally:
        if 'zip_file' in response:
            # Removing residual files
            utils.remove_files(files=files)

    return response

# ==============================================================================================================
def get_alignment_data(
        file1:str, 
        file2:str, 
        gap_pen:int=-2, 
        match_point:int=1, 
        match_pen:int=-1, 
        ignore:bool=False, 
        file_type:str='txt'):
    '''
    Extracts the header and sequence data from the input file and counts the totals in the sequence
    
    Parameter(s):
        file1 (str): path to the file containing the first sequence data
        file2 (str): path to the file containing the second sequence data
        gap_pen (int, default=-2): penalty for gaps in the alignment
        match_point (int, default=1): point(s) added to the score for matches in the alignment
        match_pen (int, default=-1): penalty for mismatches in the alignment
        ignore (bool, default=False): condition to ignore start and end gap penalties for local alignment
        file_type (str, default=txt): type of file(s) to be returned in the zip file
    
    Output(s):
        Creates a zip file composed of a file depending on file_type which contains the alignment results of a
        sequence if no errors occur. Returns a dictionary response containing the path of the zip file or any errors.
    '''

    response = {}
    files = []

    try:
        seq1 = utils.get_data(os.path.join(PATH, file1))
        seq2 = utils.get_data(os.path.join(PATH, file2))

        data = al.Alignment(
            seq=seq1[0][0], 
            ref=seq2[0][0], 
            gap_pen=gap_pen, 
            match_point=match_point, 
            match_pen=match_pen, 
            ignore=ignore
        )

        # Create a text file for each alignment type
        if file_type == 'txt':
            filename = f'./temp/alignment.txt'
            files.append(utils.make_txt(data=data.results, filename=filename))

        # Create a csv file for each alignment type
        elif file_type == 'csv':
            filename = f'./temp/alignment.csv'
            files.append(utils.make_csv(data=data.results, filename=filename))

        # Create a json file for the alignment data
        elif file_type == 'json':
            filename = f'./temp/alignment.json'
            files.append(utils.make_json(data=data.results, filename=filename))

        # Zipping the collection of files
        response['zip_file'] = utils.create_zip(files=files, zipname='./temp/alignment.zip')

    except Exception as e:
        response['error'] = str(e)
    finally:
        if 'zip_file' in response:
            utils.remove_files(files=files)

    return response

# ==============================================================================================================
def get_variance_data(
        file:str, 
        plot:bool, 
        smooth:bool, 
        vRegion:bool, 
        spacing:int=30, 
        numV:int=6,
        file_type:str='txt'):
    '''
    Calculates the variance between the sequences and plots the data
    
    Parameter(s):
        file (str): path to the file containing the sequence data
        plot (bool): save the raw/smooth plot to a file if true
        smooth (bool): smooth the data if true, else plot the raw data
        vRegion (bool): save the plot with variance regions to a file if true
        spacing (int, default=30): specifies the minimum number of bases needed for a variance region (peak)
        numV (int, default=6): specifies the minimum number of desired variance regions (peaks)
        file_type (str, default=txt): type of file(s) to be returned in the zip file
    
    Output(s):
        zipPth (str): a path to a zip file containing the ploted data
    '''

    response = {}
    files = []

    try:
        # These are hard coded (user can't change)
        file1 = './temp/plot.jpg'
        file2 = './temp/v_regions_plot.jpg'

        seq = utils.get_data(os.path.join(PATH, file))
        data = vr.Variance(sequences=seq[0], header=seq[1])

        # Plot smooth/raw data and save to a file
        if plot or smooth:
            data.plot_data(
                display=False,
                smooth=smooth,
                filename=file1
            )

            files.append(file1)

        # Plot variance regions and save to a file
        if vRegion:
            data.plot_v_regions(
                display=False,
                spacing=spacing,
                numV=numV,
                filename=file2
            )

            files.append(file2)

        # Create zip file for export
        response['zip_file'] = utils.create_zip(files=files, zipname='./temp/variance.zip')

    except Exception as e:
        response['error'] = str(e)
    finally:
        if 'zip_file' in response:
            # Delete individual files that are no longer needed
            utils.remove_files(files=files)

    return response

# ==============================================================================================================
def main():
    #sequence_data = get_sequence_data(file='testing.fna', codon='codon', amino='amino')
    variance_data = get_variance_data(file='sequences.fna', plot=True, smooth=True, vRegion=True)


if __name__ == "__main__":
    main()
