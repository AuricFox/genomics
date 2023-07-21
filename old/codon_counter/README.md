## CSCI 5481,  Homework 1

Samuel Kitzerow, kitze012

### Running Program

```
python .\count_codons.py input_file output_file option
python .\count_codons.py input_file output_file codons
python .\count_codons.py input_file output_file amino
```

Codon files:  
```
python .\count_codons.py SARS-CoV-2_separate_genes.fna ./csv_files/SARS-CoV-2_separate_genes.csv codons
python .\count_codons.py SARS-CoV-2_whole_genome.fna ./csv_files/SARS-CoV-2_whole_genome.csv codons
```

Amino acid files:  
```
python .\count_codons.py SARS-CoV-2_separate_genes.fna ./csv_files/separate_amino_acids.csv amino
python .\count_codons.py SARS-CoV-2_whole_genome.fna ./csv_files/whole_amino_acids.csv amino
```

Test files:
```
python .\count_codons.py test.fna ./csv_files/test.csv codons
python .\count_codons.py test.fna ./csv_files/test_amino.csv amino
```

Merging Files:
```
python .\setup.py merge separate_file whole_file combined_file
python .\setup.py merge ./csv_files/SARS-CoV-2_separate_genes.csv ./csv_files/SARS-CoV-2_whole_genome.csv ./csv_files/combined_codons.csv
python .\setup.py merge ./csv_files/separate_amino_acids.csv ./csv_files/whole_amino_acids.csv ./csv_files/combined_amino_acids.csv
```

### Other Files

The program `codon_mapping.py` creates an object that takes in the codon data and counts the codons and amino acids. The program `setup.py` was used just to setup the matrix `self.amino` for mapping codons to their respective amino acids. It is also used to merge the corresponding csv files for easier processing in an excel spread sheet. It is not used in the main program.

### Processing Data

The header and empty lines are ignored when the file is read in. The line containing the genomic sequence is then passed to a `map_codons` object. Three characters are accessed each iteration from front to end. The three characters (codon) are then compared to corresponding ingeter values (A = 0, C = 1, G = 2, T = 3). The values are used as index locations in the matrix `self.matrix`. For example, the codon ACT evaluates to 013 or matrix 0, column 1, and row 3. In that location is a index value for the matching codon in the list `self.codon_count` which stores the number of codons. Once matched, the count for that codon will be incremented by 1. This can be done with one codon or a whole sequence. Any remaining characters are ignored

There is also an option to compute the count for amino acids. An equal length list `self.amatrix` stores the index locations for the codon's corresponding amino acids in `self.amino`. For example, the codon ACT will return the value 2 which is the index location for the amino acid Thr. The count for the amino acid is then updated in `amino_count` which is returned when all amino acids have been counted. 