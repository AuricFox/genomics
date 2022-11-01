## CSCI 5481

Samuel Kitzerow, kitze012

### Running Program

Counting codons and writing to files:  
```
python .\main.py -c input_file output_file
python .\main.py -c ./fna_files/SARS-CoV-2_separate_genes.fna ./output_files/SARS-CoV-2_separate_genes.csv codons
python .\main.py -c ./fna_files/SARS-CoV-2_whole_genome.fna ./output_files/SARS-CoV-2_whole_genome.csv codons
```

Counting amino acids and writing to files:  
```
python .\main.py -a input_file output_file
python .\main.py -a ./fna_files/SARS-CoV-2_separate_genes.fna ./output_files/separate_amino_acids.csv amino
python .\main.py -a ./fna_files/SARS-CoV-2_whole_genome.fna ./output_files/whole_amino_acids.csv amino
```

Aligning sequences with defualt penalties at start and end bases:
```
python .\main.py -l reference_fna sequence_fna output_file
python .\main.py -l ./fna_files/sars_spike_protein.fna ./fna_files/moderna_mrna.fna ./output_files/moderna_align.txt
python .\main.py -l ./fna_files/sars_spike_protein.fna ./fna_files/pfizer_mrna.fna ./output_files/pfizer_align.txt
```

Alingment of amino acid sequences with default penalties:
```
python .\main.py -la fna1 fna2 [start/end penalties] [gap penalty] [mismatch penalty] output
python .\main.py -la ./fna_files/sars_spike_protein.fna ./fna_files/moderna_mrna.fna f -2 -1 ./output_files/moderna_aa_align.txt
python .\main.py -la ./fna_files/sars_spike_protein.fna ./fna_files/pfizer_mrna.fna f -2 -1 ./output_files/pfizer_aa_align.txt
```

Aligning sequences with defualt penalties at start and end bases:
```
python .\main.py -lo fna1 fna2 [start/end penalties] [gap penalty] [mismatch penalty] output
python .\main.py -lo ./fna_files/sars_spike_protein.fna ./fna_files/moderna_mrna.fna t -2 -1 ./output_files/moderna_es_align.txt
python .\main.py -lo ./fna_files/sars_spike_protein.fna ./fna_files/pfizer_mrna.fna t -2 -1 ./output_files/pfizer_es_align.txt
```

Test alignments without writing to files:
```
python program -t fna1 fna2 [start/end penalties] [gap penalty] [mismatch penalty]
python .\main.py -t ./fna_files/sars_spike_protein.fna ./fna_files/moderna_mrna.fna t -2 -1
python .\main.py -t ./fna_files/sars_spike_protein.fna ./fna_files/pfizer_mrna.fna f -2 -1
```

Merging Files:
```
python .\setup.py merge separate_file whole_file combined_file
python .\setup.py merge ./csv_files/SARS-CoV-2_separate_genes.csv ./csv_files/SARS-CoV-2_whole_genome.csv ./csv_files/combined_codons.csv
python .\setup.py merge ./csv_files/separate_amino_acids.csv ./csv_files/whole_amino_acids.csv ./csv_files/combined_amino_acids.csv
```

### Other Files

The program `sequence.py` creates an object that takes in the codon data and counts the codons and amino acids. The program `setup.py` was used just to setup the matrix `self.amino` for mapping codons to their respective amino acids. It is also used to merge the corresponding csv files for easier processing in an excel spread sheet. It is not used in the main program.

### Processing Data

The header and empty lines are ignored when the file is read in for simple codon counting. It is used for sequence alignment. The line containing the genomic sequence is then passed to a `sequence` object. Three characters are read in at a time and compared to the codons in the dictionary `self.codon`. They are then added to the count. Any remaining characters are ignored. The same process applies to amino acids.

Once the sequence data is collected, it can then be passed to an alignment object in `alignment.py`. This program has the ability to alignment two sequences using the Needleman–Wunsch algorithm. Start and end gaps can be ignored by zeroing the initial columns and rows.