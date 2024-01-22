# Bioinformatics

Bioinformatics is an interdisciplinary field that combines biology, computer science, and mathematics to analyze and interpret biological data, particularly from large-scale genomic and molecular studies. It involves the use of computational tools and techniques to understand biological processes, genetic variations, and the relationships between genes, proteins, and other biomolecules.

## Features

### Sequence (Codon Counter)

Analizes the number of codons (3-mers), amino acids, or specified k-mers within a genetic sequence. Takes in a fna or fastq file containing genetic sequences and analizes the data within the file. The user can specify which type they want counted and/or the k number for k-mers (3-mer, 4-mer, etc.). Furthermore, the user can request certain file types to be returned (txt, csv, or json).

File(s) Returned:  
```
codons_count.(txt or csv)
amino_count.(txt or csv)
kmer_count.(txt or csv)
sequenece_count.json
```

### Sequence Alignment

Performs an allignment of two sequences to see how they compare. Takes in two fna or fastq files containing genetic sequences and aligns them using the Needleman-Wunsch algorithm. The user can specify whether they want global alignment (count end gaps) or loacl alignment (ignore end gaps). The user can also request certain file types to be returned (txt, csv, or json).

File(s) Returned:  
```
alignment.(txt, csv, json)
```

### Sequence Variance

Identifies variable regions in amplicon sequences, pieces of DNA/RNA that is the source of amplification or replication events. Takes in a fna or fastq file containing genetic sequences and plots the variance regions. The user can also request certain image file types to be returned (pdf, png, or jpg).

File(s) Returned:  
```
plot.(pdf, png, or jpeg)
v_regions_plot.(pdf, png, or jpeg)
```

### Sequence Assembly

Assembles sequence fragments into a complete genome. Takes in a fna or fastq file containing genetic fragments (reads) and assembles them using a de Bruijn graph and eulerian cycle algorithm. Users can specify the k-mer size and number of reads to assemble at a time.

Process Frow:  
De Bruijn Garph (1) -> Eulerian Cycle/Path (2) -> Alignment (3)

(1) reads.FASTQ -> Graph.pdf, edges.txt, directed_graph.txt
        Converts reads into kmers and edges. This data is then used to constuct a directed graph file and figure.

(2) directed_graph.txt -> eulerianPath.txt, eulerianCycle.txt (assembled contigs)
        Converts a directed graph file into an Eulerian cycle/path file. This is an assembled contig.

(3) eulerianPath.txt, eulerianCycle.txt -> align.txt, comparison.pdf
        Compares the Eulerian cycle/path file with the assembled spike protein and creates text file with the comparisons and a plot.

File(s) Returned:  
```
de_bruijn_graph.pdf
directed_edges.txt
edges.txt
alignment_plot.pdf
sequence_alignment.txt
```

### Phylogeny

Construct a phylogeny tree from a series of genetic sequences. Takes in a fna or fastq file containing genetic sequences and implements the Nei-Saitou neighbor-joining algorithm for phylogeny construction. The final results are then returned.

File(s) Returned:  
```
edges.txt
genetic-distances.txt
tree.pdf
```

## Getting Started

To get started with Genomics, follow these steps:

1. **Clone the Repository:**
    ```
    git clone https://github.com/AuricFox/genomics.git
    ```

2. **Navigate to the Project Directory:**
    ```
    cd genomics
    ```

3. **Setup Environment:**
    ```
    pip install virtualenv  
    virtualenv env

    .\env\Scripts\activate      # Windows
    source env/bin/activate     # Mac OS
    ```

4. **Install Dependencies:**
    ```
    (env) pip install flask
    (env) pip install numpy
    (env) pip install networkx
    (env) pip install matplotlib
    (env) pip install typing
    (env) pip install Bio
    (env) pip install logging
    (env) pip install minetypes
    ```

5. **Run Server:**
    ```
    (env) python app.py
    ```

    The server will start running, and you can access the application by navigating to `http://localhost:5000` in your web browser.
