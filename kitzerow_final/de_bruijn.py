# De Bruijn Graph
'''
Sources
https://www.nature.com/articles/nbt.2023
https://eaton-lab.org/slides/genomics/answers/nb-10.2-de-Bruijn.html
'''

import toyplot as ty

class De_bruijn:
    def __init__(self, seq = [], header = []):
        self.header = header
        self.seq = seq

        self.kmers = {}
        self.edges = set()
        self.de_bruijn_graph()

    # ----------------------------------------------------------------------------------------------------------
    # Sets values for attributes self.kmer and self.edges
    def de_bruijn_graph(self, k=3, cycle=True):

        self.kmers = self.get_kmers(k, cycle)
        self.edges = self.get_edges(self.kmers)

    # ----------------------------------------------------------------------------------------------------------
    # Build a list of all kmers in the provided sequences
    # k: kmer size
    # cycle: sequence is cyclic or not
    def get_kmers(self, k=3, cycle=True):
        kmers = {}

        for s in self.seq:                          # Iterate thru all sequences
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
        
                if kmer in kmers:              # Add to kmer count (kmer is already found)
                    kmers[kmer] += 1
                else:                               # Add kmer to dictionary (kmer is new)
                    kmers[kmer] = 1
        
        return kmers

    # ----------------------------------------------------------------------------------------------------------
    # Finds the edges of the kmers
    def get_edges(self, kmers):
        edges = set()

        for k1 in kmers:
            for k2 in kmers:
                if k1 != k2:                                # Iterate thru all non-equal kmers (k1 != k2)
                    
                    if k1[1:] == k2[:-1]:                   # k-1 mers are the same (ex. k1=ACGT, k2=CGTA, CGT == CGT)
                        edges.add((k1[:-1], k2[:-1]))       # Add edge

                    if k1[:-1] == k2[1:]:                   # k-1 mers are the same (ex. k1=CGTA, k2=ACGT, CGT == CGT)
                        edges.add((k2[:-1], k1[:-1]))       # Add edge

        return edges

    # ----------------------------------------------------------------------------------------------------------
    # Creates graph of connected nodes with corresponding kmers
    def plot_graph(self, width=500, height=500):
    
        graph = ty.graph(
            [i[0] for i in self.edges],
            [i[1] for i in self.edges],
            width=width,
            height=height,
            tmarker=">", 
            vsize=25,
            vstyle={"stroke": "black", "stroke-width": 2, "fill": "none"},
            vlstyle={"font-size": "11px"},
            estyle={"stroke": "black", "stroke-width": 2},
            layout=ty.layout.FruchtermanReingold(edges=ty.layout.CurvedEdges()))
        return graph
    