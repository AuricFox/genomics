# De Bruijn Graph
'''
Sources
https://www.nature.com/articles/nbt.2023
https://eaton-lab.org/slides/genomics/answers/nb-10.2-de-Bruijn.html
'''

import matplotlib.pyplot as plt
import networkx as nx               # Plotting directed graph
from tqdm import tqdm               # Progress Bar
import math

class De_bruijn:
    def __init__(self, seq = [], header = []):
        self.header = header
        self.seq = seq

        self.kmers = {}
        self.edges = set()

    # ----------------------------------------------------------------------------------------------------------
    def de_bruijn_graph(self, start=0, end=1, k=3, cycle=True):

        if(start < 0 or start > end or start == end):   # Cannot exceed bounds of the list self.seq
            print("ERROR: Invalid input!")
            return

        if(end > len(self.seq)):                        # End exceeds sequence list, change it to list length
            end = len(self.seq)

        seq = self.seq[start:end]                       # Get kmers from this list of sequences
        self.kmers = self.get_kmers(seq, k, cycle)
        self.edges = self.get_edges(self.kmers)

    # ----------------------------------------------------------------------------------------------------------
    # Build a list of all kmers in the provided sequences
    # k: kmer size
    # cycle: sequence is cyclic or not
    def get_kmers(self, seq, k=3, cycle=False):
        kmers = {}

        for s in tqdm(seq, desc=str(k)+'-mers'):   # Iterate thru sequences with progress bar
            #print("Sequence: ", s)

            for i in range(0, len(s)):              # Iterate thru a sequence to find kmers
                kmer = s[i:i + k]

                length = len(kmer)
                if cycle:                           # Get kmer for cyclic sequences
                    if len(kmer) != k:
                        kmer += s[:(k - length)]    # Cyclic kmer (ex. ACCGATCG, cyclic 3-mer = CGA and GAC)

                else:                               # Skip for non-cyclic sequences
                    if len(kmer) != k:
                        continue

                if kmer in kmers:                   # Add to kmer count (kmer is already found)
                    kmers[kmer] += 1
                else:                               # Add kmer to dictionary (kmer is new)
                    kmers[kmer] = 1
        
        return kmers

    # ----------------------------------------------------------------------------------------------------------
    # Finds the edges of the kmers
    def get_edges(self, kmers):
        edges = set()

        for k1 in tqdm(kmers, desc='Edges'):
            for k2 in kmers:
                if k1 != k2:                                # Iterate thru all non-equal kmers (k1 != k2)

                    if k1[1:] == k2[:-1]:                   # k-1 mers are the same (ex. k1=ACGT, k2=CGTA, CGT == CGT)
                        edges.add((k1[:-1], k2[:-1]))       # Add edge

                    if k1[:-1] == k2[1:]:                   # k-1 mers are the same (ex. k1=CGTA, k2=ACGT, CGT == CGT)
                        edges.add((k2[:-1], k1[:-1]))       # Add edge

        return edges

    # ----------------------------------------------------------------------------------------------------------
    # Create edges.txt file for creation of spike_protein_directed_graph.txt
    def create_edges_file(self, file='output/temp/edges.txt'):
        print("Creating Edge File: ", file)

        with open(file, 'w') as f:
            for edge in tqdm(self.edges, desc='Edge File'):
                x1, x2 = edge
                f.write(x1 + '->' + x2 + '\n')
        f.close()

    # ----------------------------------------------------------------------------------------------------------
    # Writes edge data to a text file
    def create_directed_graph(self, file='output/temp/spike_protein_directed_graph.txt'):
        print("Creating Directed Graph File: ", file)

        with open(file, "w") as f:
            added_nodes = set()                                     # Set of (nodes, destination) pairs that have already been iterated through

            for edge in tqdm(self.edges, desc='Directed Graph'):
                node,dest = edge
                if edge not in added_nodes:                         # Check if (node, destination) pair has already been writted
                    f.write(node + ' -> ' + dest)

                    for edge2 in self.edges:                        # Iterate through all edges to see if there are any other edges coming out of node
                        node2, dest2 = edge2
                        # If the nodes are the same AND (node, destination) pair is not the same as above AND this node has not already been created
                        if node == node2 and edge != edge2 and edge2 not in added_nodes:
                            f.write(',' + dest2)                    # Write the additional destination for exising node
                            added_nodes.add(edge2)                  # add (node, destination) pair to already added set

                    f.write('\n')
        f.close()

    # ----------------------------------------------------------------------------------------------------------
    # Creates graph of connected nodes with corresponding kmers using matplotlib
    def matplot_graph(self, show_lab=False, show_fig=True, save_fig=False, file='./output/graph/deBruijn.png'):
        print("Creating Garph Image: ", file)
        with tqdm(total=4, desc='Image Graph') as bar:              # Progress bar
            plt.clf()
            fig = nx.DiGraph()                  # Initialize weighted graph
            fig.add_edges_from(self.edges)      # Add edges to weighted graph
            bar.update(1)
            k = 0.9/math.sqrt(fig.order())
            pos = nx.shell_layout(fig, scale=2)          # Set layout to shell (circular)
            bar.update(1)

            options = {
                "node_color": "#A0CBE2",
                "node_size": 20,
                "edge_color": "#7d0901",
                "with_labels": show_lab,            # Show kmer labels
                "font_size": 5,
                "font_color": "#0a0a0a"
            }

            nx.draw(fig, pos, **options)        # Drawing directed graph
            bar.update(1)
            plt.axis("off")                     # Do not show any axis
            if(show_fig): plt.show()            # Toggel show
            if(save_fig): plt.savefig(file, transparent=False ,dpi=500)     # Saving file and setting size
            bar.update(1)

    # ----------------------------------------------------------------------------------------------------------
    # Master function for creating documents
    # edge_graph: creates graph image
    # edge_file: creates edge file used for directed graph
    # dir_graph: creates directed edge file
    # i: label for k-mer used (used for kmer loop)
    def make_docs(self, edge_graph=False, edge_file=False, dir_graph=False, i='1'):

        if(edge_graph): self.matplot_graph(False, False,True, './output/graph/deBruijn_' + i + '.png')
        if(edge_file): self.create_edges_file('output/temp/edges_' + i + '.txt')

        if(dir_graph): 
            file = 'output/temp/spike_protein_directed_graph_' + i + '.txt'
            self.create_directed_graph(file)
            return file
