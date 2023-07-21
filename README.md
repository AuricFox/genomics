# Genomics
Repo dedicated to a computer science class on genomics.

## Sequence (Codon Counter)

The program `codon_mapping.py` creates an object that takes in the codon data and counts the codons and amino acids.

## Sequence Alignment

## Sequence Variance

## Sequence Assembly

This python script assembles genetic sequences from read fragments. The de Bruijn graph breaks up the reads into
individual codons and edges to form a graph. Then a eulerian cycle is found.

Process Frow:
De Bruijn Garph -> Eulerian Cycle/Path -> Alignment  
     (1)        ->         (2)         ->    (3)

(1) reads.FASTQ -> Graph.pdf, edges.txt, directed_graph.txt
        Converts reads into kmers and edges. This data is then used to constuct a directed graph file and figure.

(2) directed_graph.txt -> eulerianPath.txt, eulerianCycle.txt (assembled contigs)
        Converts a directed graph file into an Eulerian cycle/path file. This is an assembled contig.

(3) eulerianPath.txt, eulerianCycle.txt -> align.txt, comparison.pdf
        Compares the Eulerian cycle/path file with the assembled spike protein and creates text file with the comparisons and a plot.

File(s):
./input/sars_spike_protein_reads.fastq

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