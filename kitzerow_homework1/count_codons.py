# Samuel Kitzerow, kitze012
# Homework 1, Parsing fna file and making a csv file
import csv
import codon_mapping as map

# Wrights the data to a csv file
def make_csv(filename, data):
    f = filename.split(".")     # Using the same filename but as a csv
    filename = f[0] + ".csv"
    print(filename, data)
    
    file = open()
    with open(filename, 'w') as file:
        cfile = csv.writer(file)
        cfile.writerow(['Codon', 'Count'])  # Wrights header to csv
        cfile.writerows(data)               # Wrights codon data to csv

# =======================================================================================
# Retreives data from the fna file and parses it
def get_data(filename):
    file = open(filename, 'r')

    genes = map.map_codons()
    
    for line in file:
        line = line.strip()
        
        if(line == ""):             # Empty line at the top of the file
            # print("Not text")
            continue
        elif(line[0] == ">"):       # Header information
            # print("Header")
            continue
        else:                       # Genetic sequence containing desired codons
            # print("Codons")
            genes.add_sequence(line)
    
    print(genes.get_codon_count())

# =======================================================================================

if __name__ == "__main__":

    # get_data('SARS-CoV-2_whole_genome.fna')
    get_data('SARS-CoV-2_separate_genes.fna')

    #test = map.map_codons()
    #test.debugging()
