# Samuel Kitzerow, kitze012
# Homework 1, Parsing fna file and making a csv file
import csv
import codon_mapping as map

# Wrights the data to a csv file
def make_csv(filename, data):
    f = filename.split(".")     # Using the same filename but as a csv
    filename = f[0] + ".csv"
    print(filename, data)

    with open(filename, 'w') as file:
        cfile = csv.writer(file)
        cfile.writerow(['Codon', 'Count'])  # Wrights header to csv
        cfile.writerows(data)               # Wrights codon data to csv

# =======================================================================================
# Retreives data from the fna file and parses it
def get_data(filename):
    fd = open(filename, 'r')
    contents = fd.read()

# =======================================================================================

if __name__ == "__main__":

    # get_data('SARS-CoV-2_whole_genome.fna')
    # get_data('SARS-CoV-2_separate_genes')

    test = map.map_codons()
    test.debugging()
