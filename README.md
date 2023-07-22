# Bioinformatics

Bioinformatics is an interdisciplinary field that combines biology, computer science, and mathematics to analyze and interpret biological data, particularly from large-scale genomic and molecular studies. It involves the use of computational tools and techniques to understand biological processes, genetic variations, and the relationships between genes, proteins, and other biomolecules.

## Sequence (Codon Counter)

Analizes the number of codons (3-mers), amino acids, or specified k-mers within a genetic sequence. Takes in a fna or fastq file containing genetic sequences and analizes the data within the file. Users can specify which type they want counted and/or the k number for k-mers (3-mer, 4-mer, etc.).

## Sequence Alignment

Performs an allignment of two sequences to see how they compare. Takes in two fna or fastq files containing genetic sequences and aligns them using the Needleman-Wunsch algorithm. Users can specify whether they want global alignment (count end gaps) or loacl alignment (ignore end gaps).

## Sequence Variance

Identifies variable regions in amplicon sequences, pieces of DNA/RNA that is the source of amplification or replication events. Takes in a fna or fastq file containing genetic sequences and plots the variance regions.

## Sequence Assembly

Assembles sequence fragments into a complete genome. Takes in a fna or fastq file containing genetic fragments (reads) and assembles them using a de Bruijn graph and eulerian cycle algorithm. Users can specify the k-mer size and number of reads to assemble at a time.

Process Frow:  
De Bruijn Garph (1) -> Eulerian Cycle/Path (2) -> Alignment (3)

(1) reads.FASTQ -> Graph.pdf, edges.txt, directed_graph.txt
        Converts reads into kmers and edges. This data is then used to constuct a directed graph file and figure.

(2) directed_graph.txt -> eulerianPath.txt, eulerianCycle.txt (assembled contigs)
        Converts a directed graph file into an Eulerian cycle/path file. This is an assembled contig.

(3) eulerianPath.txt, eulerianCycle.txt -> align.txt, comparison.pdf
        Compares the Eulerian cycle/path file with the assembled spike protein and creates text file with the comparisons and a plot.

File(s): `./input/sars_spike_protein_reads.fastq`

## Phylogeny

# Server-side
## Server Setup

This server is run using the FLASK framework used in Python, But an envirnment must first be setup.

STEP 1: cd into working directory that contains your project  

STEP 2: Install env module: `pip install virtualenv`  

STEP 3: Activate env:  
```
C: virtualenv env               # env is the environment file name
C: env\Scripts\activate     # Windows
C: source env/bin/activate      # Mac
```  

STEP 4: Install flask in env:  
```
(env) pip install flask
```  

## Running Server

This is a development server so everything will be running on localhost (possibly on local network).

The env needs to be activated to run:  
```
C: env\Scripts\activate     # Windows
C: source env/bin/activate      # Mac
```

Execute server program:  
```
(env) python server.py
```

NOTE: Restart the server if any changes are made to any of the files.

Deactivate Environment:  
```
(env) deactivate
```