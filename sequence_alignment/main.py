# Samuel Kitzerow, kitze012
# Homework 1, Parsing fna file and making a csv file
import csv
import sys
import sequence as seq
import alignment as alg

# =======================================================================================
# Wrights the amino acid sequence data to a aa file
def make_aa(filename, data):

    file = filename.split("/")[-1]                  # ./fna_files/file.fna -> file.fna
    file = file.split(".")[0]                       # file.fna -> file
    file = "./output_files/" + file + ".aa"         # file -> ./output_files/file.aa

    with open(file, 'w', newline='') as f:
        for x in data:
            f.write(x + '\n')

# =======================================================================================
# Wrights the sequence alignment data to a text file
def make_txt(filename, data):

    with open(filename, 'w', newline='') as file:
        for x in data:
            file.write(x + '\n')

# =======================================================================================
# Wrights the codon and amino acid data to a csv file
def make_csv(filename, data, flag):
    # print(filename, data)
    
    with open(filename, 'w', newline='') as file:
        cfile = csv.writer(file)

        if (flag == 'codon'):
            cfile.writerow(['Codon', 'Count'])          # Wrights codon header to csv
        elif (flag == 'amino_acid'):
            cfile.writerow(['Amino Acid', 'Count'])     # Wrights amino acid header to csv
        else:
            print("ERROR: INVALID TYPE, codon/amino_acid ONLY!")
            return

        cfile.writerows(data)                           # Wrights codon or amino acid data to csv

# =======================================================================================
# Retreives data from the fna file and parses it
def get_data(filename, header=False):
    file = open(filename, 'r')
    seq_data = seq.sequence()
    
    for line in file:                                   # Read each line in file
        line = line.strip()                             # Strip newline characters
        
        if(line == ""):                                 # Empty line at the top of the file
            # print("Not text")
            continue
        elif(line[0] == ">" and header == True):        # Enter header information
            seq_data.add_header(line)
        elif(line[0] == ">" and header == False):       # No header information
            continue
        else:                                           # Genetic sequence containing desired codons
            seq_data.add_to_sequence(line)              # Add sequence to genes
    
    file.close()
    return seq_data

# =======================================================================================

def main():
    # ----------------------------------------------------------------------------------------------------------
    # Processes codon data from fna file and writes codons with their respective counts to a csv file
    # python .\main.py -c input_file output_file
    if (len(sys.argv) == 4 and sys.argv[1] == "-c"):
        fna_file = sys.argv[2]                                  # Extracting fna file name from arguments
        csv_file = sys.argv[3]                                  # Extracting csv file name from arguments

        genes = get_data(fna_file)                              # Processing genome data from fna file
        make_csv(csv_file, genes.get_codon_count(), 'codon')    # Writing genome data to csv file
    # ----------------------------------------------------------------------------------------------------------
    # Processes codon data from fna file and writes amino acids with their respective counts to a csv file
    # python .\main.py -a input_file output_file
    elif(len(sys.argv) == 4 and sys.argv[1] == "-a"):                       # Getting amino acid details
        fna_file = sys.argv[2]
        csv_file = sys.argv[3]

        amino_acids = get_data(fna_file)
        make_csv(csv_file, amino_acids.get_amino_count(), 'amino_acid')     # Writing amino acid data to csv file
    # ----------------------------------------------------------------------------------------------------------
    # Processes two sequences for alignment and prints the results to a text file (w/ default start and end penalties)
    # python .\main.py -l reference_fna sequence_fna output_file
    elif(len(sys.argv) == 5 and sys.argv[1] == "-l"):
        refernce_fna = sys.argv[2]
        sequence_fna = sys.argv[3]
        txt_file = sys.argv[4]

        ref_seq = get_data(refernce_fna, True)                      # Setting reference sequence
        align_seq = get_data(sequence_fna, True)                    # Setting sequence 1
        a = alg.alignment(ref_seq.sequence, align_seq.sequence)     # Setting alignment of two sequences
        a_seq = a.get_alignment()                                   # Getting alignment (ref_seq, vis, align_seq, score)

        print("Match: ",a.counts[0], " Mismatch: ",a.counts[1], " Gaps: ", a.counts[2], " E/S Gaps: ", a.counts[3])
        data = [str(a_seq[3]), ref_seq.header, a_seq[0], a_seq[1], a_seq[2], align_seq.header]
        make_txt(txt_file, data)
    
    # ----------------------------------------------------------------------------------------------------------
    # Processes two amino acid sequences for alignment and prints the results to a text file and aa file
    # python program -la fna1 fna2 [start/end penalties] [gap penalty] [mismatch penalty] output
    elif(len(sys.argv) == 8 and sys.argv[1] == "-la"):
        refernce_fna = sys.argv[2]
        sequence_fna = sys.argv[3]
        gap_pen = int(sys.argv[5])
        match_pen = int(sys.argv[6])
        txt_file = sys.argv[7]

        if(sys.argv[4] == "t"): ign = True                                      # Ignore start and end gaps
        elif(sys.argv[4] == "f"): ign = False                                   # Track start and end gaps
        else:
            print("ERROR: INVALID INPUT, t/f ONLY!")
            return

        ref_seq = get_data(refernce_fna, True)                                  # Setting reference sequence
        align_seq = get_data(sequence_fna, True)                                # Setting sequence 1

        make_aa(refernce_fna, [ref_seq.header, ref_seq.codon_to_amino()])       # Creating aa file for reference sequence
        make_aa(sequence_fna, [align_seq.header, align_seq.codon_to_amino()])   # Creating aa file for sequence 1

        a = alg.alignment(ref_seq.codon_to_amino(), align_seq.codon_to_amino(), gap_pen, match_pen, ign)    # Setting alignment of two sequences w/penalties
        a_seq = a.get_alignment()                                                                           # Getting alignment (ref_seq, vis, align_seq, score)

        print("Match: ",a.counts[0], " Mismatch: ",a.counts[1], " Gaps: ", a.counts[2], " E/S Gaps: ", a.counts[3])
        data = [str(a_seq[3]), ref_seq.header, a_seq[0], a_seq[1], a_seq[2], align_seq.header]
        make_txt(txt_file, data)

    # ----------------------------------------------------------------------------------------------------------
    # Processes two sequences for alignment with penalties and prints the results to a text file
    # python program -lo fna1 fna2 [start/end penalties] [gap penalty] [mismatch penalty] output
    elif(len(sys.argv) == 8 and sys.argv[1] == "-lo"):
        refernce_fna = sys.argv[2]
        sequence_fna = sys.argv[3]
        gap_pen = int(sys.argv[5])
        match_pen = int(sys.argv[6])
        txt_file = sys.argv[7]

        if(sys.argv[4] == "t"): ign = True                          # Ignore start and end gaps
        elif(sys.argv[4] == "f"): ign = False                       # Track start and end gaps
        else:
            print("ERROR: INVALID INPUT, t/f ONLY!")
            return

        ref_seq = get_data(refernce_fna, True)                                              # Setting reference sequence
        align_seq = get_data(sequence_fna, True)                                            # Setting sequence 1
        a = alg.alignment(ref_seq.sequence, align_seq.sequence, gap_pen, match_pen, ign)    # Setting alignment of two sequences w/penalties
        a_seq = a.get_alignment()                                                           # Getting alignment (ref_seq, vis, align_seq, score)

        print("Match: ",a.counts[0], " Mismatch: ",a.counts[1], " Gaps: ", a.counts[2], " E/S Gaps: ", a.counts[3])
        data = [str(a_seq[3]), ref_seq.header, a_seq[0], a_seq[1], a_seq[2], align_seq.header]
        make_txt(txt_file, data)
    
    # ----------------------------------------------------------------------------------------------------------
    # Run test prints but don't make csv or text files
    # python program -t fna1 fna2 [start/end penalties] [gap penalty] [mismatch penalty]
    elif(len(sys.argv) == 7 and sys.argv[1] == "-t"):
        refernce_fna = sys.argv[2]
        sequence_fna = sys.argv[3]
        gap_pen = int(sys.argv[5])
        match_pen = int(sys.argv[6])

        if(sys.argv[4] == "t"): ign = True                          # Ignore start and end gaps
        elif(sys.argv[4] == "f"): ign = False                       # Track start and end gaps
        else:
            print("ERROR: INVALID INPUT, t/f ONLY!")
            return

        ref_seq = get_data(refernce_fna, True)                                              # Setting reference sequence
        align_seq = get_data(sequence_fna, True)                                            # Setting sequence 1
        a = alg.alignment(ref_seq.sequence, align_seq.sequence, gap_pen, match_pen, ign)    # Setting alignment of two sequences w/penalties
        a_seq = a.get_alignment()                                                           # Getting alignment (ref_seq, vis, align_seq, score)

        print("Match: ",a.counts[0], " Mismatch: ",a.counts[1], " Gaps: ", a.counts[2], " E/S Gaps: ", a.counts[3])
        data = [str(a_seq[3]), ref_seq.header, a_seq[0], a_seq[1], a_seq[2], align_seq.header]

        for x in data:
            print(x)

        # ----------------------------------------------------------------------------------------------------------

    else:
        print("ERROR: INVALID ARGUMENMTS!")

# ================================================================================================================================================
if __name__ == "__main__":
    main()
    
