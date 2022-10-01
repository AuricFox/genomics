# Samuel Kitzerow, kitze012
# Homework 1, Parsing fna file and making a csv file
import csv
import sys
import codon_mapping as map

# Wrights the data to a csv file
def make_csv(filename, data, flag):
    # print(filename, data)
    
    with open(filename, 'w', newline='') as file:
        cfile = csv.writer(file)

        if (flag == 'c'):
            cfile.writerow(['Codon', 'Count'])          # Wrights codon header to csv
        elif (flag == 'a'):
            cfile.writerow(['Amino Acid', 'Count'])     # Wrights amino acid header to csv

        cfile.writerows(data)                           # Wrights codon or amino acid data to csv

# =======================================================================================
# Retreives data from the fna file and parses it
def get_data(filename):
    file = open(filename, 'r')

    genes = map.map_codons()
    
    for line in file:                   # Read each line in file
        line = line.strip()             # Strip newline characters
        
        if(line == ""):                 # Empty line at the top of the file
            # print("Not text")
            continue
        elif(line[0] == ">"):           # Header information
            # print("Header")
            continue
        else:                           # Genetic sequence containing desired codons
            # print("Codons")
            genes.add_sequence(line)    # Add sequence to genes
    
    #print(genes.get_codon_count())
    return genes                        # Return genes object

# =======================================================================================

if __name__ == "__main__":
    '''
    part_genes = get_data('SARS-CoV-2_separate_genes.fna')
    make_csv('SARS-CoV-2_separate_genes.csv', part_genes.get_codon_count())
    
    whole_genes = get_data('SARS-CoV-2_whole_genome.fna')
    make_csv('SARS-CoV-2_whole_genome.csv', whole_genes.get_codon_count())
    '''

    if (len(sys.argv) == 4 and sys.argv[3] == "-c"):        # Getting codon detail
        fna_file = sys.argv[1]                                  # Extracting fna file name from arguments
        csv_file = sys.argv[2]                                  # Extracting csv file name from arguments
        # print(fna_file, csv_file)

        genes = get_data(fna_file)                              # Processing genome data from fna file
        make_csv(csv_file, genes.get_codon_count(), 'c')        # Writing genome data to csv file

    elif(len(sys.argv) == 4 and sys.argv[3] == "-a"):        # Getting amino acid details
        fna_file = sys.argv[1]
        csv_file = sys.argv[2]

        amino_acids = get_data(fna_file)
        make_csv(csv_file, amino_acids.get_amino_count(), 'a')  # Writing amino acid data to csv file

    else:
        print("ERROR: INVALID ARGUMENMTS!")


    
    
