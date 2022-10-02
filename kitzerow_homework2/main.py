# Samuel Kitzerow, kitze012
# Homework 1, Parsing fna file and making a csv file
import csv
import sys
import sequence as map

'''
23
>SARS-Cov-2_reference_genome_spike_protein
__C_GATT_TCGGG...
  |  ||| |xx||
AGCA_ATTCTTCGG...
>Pfizer_mRNA_vaccine
'''
# =======================================================================================
# Wrights the sequence alignment data to a text file
def sequence_to_tfile(filename, data):

    with open(filename, 'w', newline='') as file:
        file.writelines(data)

# =======================================================================================
# Wrights the codon and amino acid data to a csv file
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
def get_data(filename, flag):
    file = open(filename, 'r')
    genes = map.sequence()
    header = ""
    
    for line in file:                                   # Read each line in file
        line = line.strip()                             # Strip newline characters
        
        if(line == ""):                                 # Empty line at the top of the file
            # print("Not text")
            continue
        elif(line[0] == ">" and flag == "-l"):          # Header information
            header = line
        elif(line[0] == ">" and flag != "-l"):          # Header information
            continue
        elif(flag == "-a" or flag == "-c"):
            genes.add_to_count(line)                    # Add sequence to genes
        else:                                           # Genetic sequence containing desired codons
            # print("Codons")
            genes.add_to_sequence(line)                 # Add sequence to genes
    
    #print(genes.get_codon_count())
    file.close()
    return genes                        # Return genes object

# =======================================================================================

if __name__ == "__main__":

    if (len(sys.argv) == 4 and sys.argv[1] == "-c"):            # Getting codon detail
        fna_file = sys.argv[1]                                  # Extracting fna file name from arguments
        csv_file = sys.argv[2]                                  # Extracting csv file name from arguments

        genes = get_data(fna_file)                              # Processing genome data from fna file
        make_csv(csv_file, genes.get_codon_count(), 'c')        # Writing genome data to csv file
    # ----------------------------------------------------------------------------------------------------------

    elif(len(sys.argv) == 4 and sys.argv[1] == "-a"):           # Getting amino acid details
        fna_file = sys.argv[1]
        csv_file = sys.argv[2]

        amino_acids = get_data(fna_file)
        make_csv(csv_file, amino_acids.get_amino_count(), 'a')  # Writing amino acid data to csv file
    # ----------------------------------------------------------------------------------------------------------

    elif(len(sys.argv) == 5 and sys.argv[1] == "-l"):           # Get sequence alignment
        # python program [flag] fna1 fna2 output
        fna_file = sys.argv[1]
        csv_file = sys.argv[2]

        amino_acids = get_data(fna_file)
        make_csv(csv_file, amino_acids.get_amino_count(), 'a')  # Writing amino acid data to csv file
    # ----------------------------------------------------------------------------------------------------------

    elif(len(sys.argv) == 8 and sys.argv[1] == "-l"):           # Get sequence alignment
        # python program [flag] fna1 fna2 [flage] [gap penalty] [mismatch penalty] output

        fna_file = sys.argv[1]
        csv_file = sys.argv[2]

        amino_acids = get_data(fna_file)
        make_csv(csv_file, amino_acids.get_amino_count(), 'a')  # Writing amino acid data to csv file
    # ----------------------------------------------------------------------------------------------------------
    # Run test prints but don't make csv or text files
    elif(len(sys.argv) == 5 and sys.argv[1] == "-t"):           # Get sequence alignment
        # python program [flag] fna1 fna2 output
        fna_file = sys.argv[1]
        csv_file = sys.argv[2]

        amino_acids = get_data(fna_file)
        
    # ----------------------------------------------------------------------------------------------------------

    else:
        print("ERROR: INVALID ARGUMENMTS!")


    
    
