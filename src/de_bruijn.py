import numpy as np
import networkx as nx
import alignment as al
import matplotlib.pyplot as plt
from typing import List

import matplotlib, math, copy, random, time, os, utils
matplotlib.use('agg')

PATH = os.path.dirname(os.path.abspath(__file__))

'''
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
'''

class De_bruijn:
    def __init__(self, sequences:List[str]=[], header:List[str]=[], k:int=3):
        '''
        Initializes the de Bruijn graph
        
        Parameter(s):
            sequences (List[str]): a series of sequences being analyzed
            header (List[str]): a series of heder info for the corresponding sequence
            k (int, default=3): size of the k-mer, length of the substrings created from the main sequence

        Output(s):
            None
        '''
        self.header = header
        self.sequences = sequences
        self.k = k

        self.kmers = {}             # Stores all the possible k-mer combinations from the input sequences
        self.edges = set()          # Stores all the possible edges connecting the k-mers
        self.dir_graph = {}         # Stores all the directed edges from the k-mers
        self.contigs = []           # Stores the ordered sequence of contigs for assembly
        self.final_sequence = ''    # The final assembled sequence

    # ----------------------------------------------------------------------------------------------------------
    def de_bruijn_graph(self, start:int=0, end:int=-1, cycle:bool=False):
        '''
        Gets the k-mers and edges from a set number of sequence reads
        
        Parameter(s):
            start (int, default=0): starting index for the list of sequence reads (self.sequences)
            end (int, default=-1): ending index for the list of sequence reads (self.sequences)
            cycle (bool, default=False): cycle back to the start of the sequence when building k-mers
        
        Output(s): None
        '''

        # Cannot exceed bounds of the list self.sequences
        if(start < 0 or start > end or start == end):
            raise utils.InvalidInput(f"Invalid Input: {start} is not less than {end} or greater than 0!")

        # End exceeds sequence list, change it to list length
        if(end > len(self.sequences)):
            end = len(self.sequences)

        # Get kmers from this list of sequences
        self.get_kmers(self.sequences[start:end], cycle)
        self.get_edges()
        self.get_directed_graph()
        self.get_eulerian_cycle()
        self.get_assembled_str()

    # ----------------------------------------------------------------------------------------------------------
    def get_kmers(self, sequences:List[str], cycle=False):
        '''
        Builds a list of all kmers (sub-strings) in the provided sequences

        Parameter(s):
            sequences (str): a list of sequences being broken up into individual k-mers/sub-strings
            cycle (bool, defualt=False): sequence is cyclic or not

        Output(s): None
        '''

        for sequence in sequences:

            for i in range(0, len(sequence)):                   # Iterate thru a sequence to find kmers
                kmer = sequence[i:i + self.k]

                length = len(kmer)
                if cycle:                                       # Get kmer for cyclic sequences
                    if len(kmer) != self.k:
                        kmer += sequence[:(self.k - length)]    # Cyclic kmer (ex. ACCGATCG, cyclic 3-mer = CGA and GAC)

                else:                                           # Skip for non-cyclic sequences
                    if len(kmer) != self.k:
                        continue

                if kmer in self.kmers:                          # Add to kmer count (kmer is already found)
                    self.kmers[kmer] += 1
                else:                                           # Add kmer to dictionary (kmer is new)
                    self.kmers[kmer] = 1

    # ----------------------------------------------------------------------------------------------------------
    def get_edges(self):
        '''
        Finds all the connecting edges of the k-mers for the de Bruijn garph. Populate the kmer_dict 
        with k-1 mers as keys and lists of corresponding kmers as values. Example, kmer = ACTG, k1 = ACT, 
        and k2 = CTG. If k1 or k2 is not in the kmer dictionary, they are added.
        
        Parameter(s): None
        
        Output(s): None
        '''

        kmer_dict = {}  # Create a dictionary to store k-1 mers and their corresponding kmers

        # Populate the kmer_dict with k-1 mers and corresponding kmers
        for kmer in self.kmers:
            k1 = kmer[:-1]  # k1 is assigned the k-1 mer of the current kmer (excluding the last character).
            k2 = kmer[1:]   # k2 is assigned the k-1 mer of the current kmer (excluding the first character)

            if k1 not in kmer_dict:         # k1 is not in the dict, initialize it
                kmer_dict[k1] = []
            kmer_dict[k1].append(kmer)      # Add the kmer to k1's list

            if k2 not in kmer_dict:         # k2 is not in the dict, initialize it
                kmer_dict[k2] = []
            kmer_dict[k2].append(kmer)      # Add the kmers to k2's list

        # Create edges based on the kmer_dict
        for key, kmers in kmer_dict.items():                    # Iterate thru each item in the dict
            if len(kmers) > 1:                                  # Ignore key's with no edges
                # Add edges without repeating them
                for i in range(len(kmers)):
                    for j in range(i + 1, len(kmers)):
                        self.edges.add((kmers[i], kmers[j]))

    # ----------------------------------------------------------------------------------------------------------
    def get_directed_graph(self):
        '''
        Creates a directed graph dictionary containing a node as a key and a list of destinations as the value
        
        Parameter(s): None
        
        Output(s): None
        '''

        for edge in self.edges:
            node, dest = edge

            if node not in self.dir_graph:
                self.dir_graph[node] = []

            self.dir_graph[node].append(dest)

    # ----------------------------------------------------------------------------------------------------------
    def get_eulerian_cycle(self):
        '''
        Finds the eulerian cycle within the directed graph. Starts by calculating the total number of edges in the 
        graph and randomly picking a starting position. Then it iterates through the directed graph until all the 
        edges have been accounted for and have been added to the final sequence.

        Parameter(s): None

        Output(s):
            An ordered list of edges (str) ready to be assembled into a complete sequence
        '''

        copy_graph = copy.deepcopy(self.dir_graph)

        num_edges = 1
        values = [val for val_list in copy_graph.values() for val in val_list]
        num_edges += len(values)

        degrees = {key: [values.count(key), len(copy_graph[key]), 0] for key in copy_graph.keys()}

        cycle = []

        # Randomly select a starting node to iterate from
        current_node = random.choice([key for key in copy_graph.keys()])

        # Run until all the edges have been accounted for
        while len(self.contigs) != num_edges:

            # Current node has out going edges and is not empty
            if copy_graph[current_node] != []:
                cycle.append(current_node)
                next_possibles = copy_graph[current_node]

                # Randomly select a new node to iterate from
                new_cn = random.choice(next_possibles)
                copy_graph[current_node].remove(new_cn)
                current_node = new_cn

            # Current node has no out going edges and is empty
            elif copy_graph[current_node] == []:
                self.contigs.insert(0, current_node)
                if not cycle:
                    break
                else:
                    current_node = cycle[-1]
                    cycle.pop()

    # ----------------------------------------------------------------------------------------------------------
    def get_assembled_str(self):
        '''
        Appends the contigs together into a final sequence string

        Parameter(s): None

        Output(s): None
        '''

        self.final_sequence = self.contigs[0]                   # Start with first contig
        for contig in self.contigs[1:]:                         # Loop thru the remaining contigs
            self.final_sequence += contig[len(contig) - 1]      # Append the last letter to the final sequence

    # ----------------------------------------------------------------------------------------------------------
    def create_edges_file(self, filename:str='./temp/edges.txt'):
        '''
        Creates an edges file for testing
        
        Parameter(s):
            filename (str, default=./temp/edges.txt): file that stores the edge data

        Output(s):
            A path to the saved file containing the edge data
        '''

        file_path = os.path.join(PATH, filename)
        print("Creating Edge File: ", file_path)

        with open(file_path, 'w') as f:
            for edge in self.edges:
                x1, x2 = edge
                f.write(f"{x1}->{x2}\n")

        return file_path

    # ----------------------------------------------------------------------------------------------------------
    def create_directed_graph_file(self, filename='./temp/directed_graph.txt'):
        '''
        Creates a directed graph dictionary containing a node as a key and a list of destinations as the value.
        Writes the results to a text file. Used for testing
        
        Parameter(s): None
        
        Output(s):
            A path to the saved file containing nodes and their corresponding list of destinations.
        '''

        file_path = os.path.join(PATH, filename)
        print("Creating Directed Graph File: ", file_path)

        # Write the directed graph data to the file
        with open(file_path, "w") as f:
            for node, destinations in self.dir_graph.items():
                f.write(node + ' -> ' + ','.join(destinations) + '\n')

        return file_path

    # ----------------------------------------------------------------------------------------------------------
    def plot_graph(self, show_label:bool=False, filename:str='./temp/deBruijn.png'):
        '''
        Creates graph of connected nodes with corresponding kmers using matplotlib

        Parameter(s):
            show_label (bool, default=False): displays the labels on the plot if True else displays no labels
            filename (str, default=./temp/deBruijn.png)

        Output(s):
            A path to the saved file containing the plotted data
        '''

        file_path = os.path.join(PATH, filename)
        print("Creating Garph Image: ", file_path)
        
        plt.clf()
        fig = nx.DiGraph()                                  # Initialize weighted graph
        fig.add_edges_from(self.edges)                      # Add edges to weighted graph
        k = 0.5/math.sqrt(fig.order())                      # Used for spacing in spring layout
        #pos = nx.spring_layout(fig, k=k)
        pos = nx.shell_layout(fig, scale=2)                 # Set layout to shell (circular)

        options = {
            "node_color": "#A0CBE2",
            "node_size": 20,
            "edge_color": "#7d0901",
            "with_labels": show_label,                      # Show kmer labels
            "font_size": 5,
            "font_color": "#0a0a0a"
        }

        nx.draw(fig, pos, **options)                        # Drawing directed graph
        plt.axis("off")                                     # Do not show any axis
        plt.savefig(file_path, transparent=False ,dpi=500)  # Saving file and setting size

        return file_path

    # ----------------------------------------------------------------------------------------------------------
    def make_docs(self, edge_file:bool=False, dir_graph:bool=False, plot_type:str=None, file_type:str='.txt'):
        '''
        Master function for creating multiple documents such as an edge file, directed edge file, 
        and a plot of the constructed graph.

        Parameter(s):
            edge_graph (bool, default=False): create a plot of the graph if True else do nothing
            edge_file (bool, default=False): create an edge file if True else do nothing
            dir_graph (bool, default=False): create a directed edge file if True else do nothing
        
        Output(s):
            A list of paths to the saved files and their corresponding data (deBruijn plot, edges file, 
            and/or directed graph file)
        '''
        files = []

        if(plot_type is not None): 
            file = self.plot_graph(show_label=True, filename=f"./temp/deBruijn_{self.k}{plot_type}")
            files.append(file)
        if(edge_file): 
            file = self.create_edges_file(filename=f"output/temp/edges_{self.k}{file_type}")
            files.append(file)

        if(dir_graph): 
            file = self.create_directed_graph_file(filename=f"./temp/directed_graph_{self.k}{file_type}")
            files.append(file)

        return files

    

# ==============================================================================================================
# Aligns assembled contig with reference genome, creates text file of alignment, and plots comparison
def align_contig(filename, kmer, logging=False):
    ref_seg = utils.get_data('./input/sars_spike_protein_assembled.fna')    # Get data for reference sequence (Spike protien)

    if(logging):                                                            # Inputs are from a file, write results to a file
        data = utils.get_data(filename)                                     # Get in data from files for alignment

        a = al.alignment(ref_seg[0][0], data[0][0], -2, -1, True)           # Initialize alignment
        a.alignment_file(kmer)
        a.plot_compare(kmer)
    else:                                                                   # Inputs are a string, return string
        a = al.alignment(ref_seg[0][0], filename, -2, -1, True)             # Initialize alignment
        return a.get_alignment()

# ==============================================================================================================
# Main Functions that run the de Bruin alogrithm
# ==============================================================================================================
# Loops thru a k-mer range, builds a set of kmers and edges
# k_i, k_f: Start and end kmer size
# l_i, l_f: Start and end index of sequence reads
# res_align: Boolean value for aligning assembled sequence with the reference sequence
# res_time: Boolean value for measuring program run-time
# res_file: Boolean value for creating files
def loop_kmer(file, k_i, k_f, l_i=0, l_f=1, res_align=False, res_time=False, logging=False):
    data = utils.get_data(file)                                     # Processing genome data from fna file
    runtime = []
    num_reads = l_f - l_i                                           # Number of reads being assembled

    db_data = De_bruijn(data[0], data[1])                           # Create de Bruijn graph

    for i in range(k_i, k_f+1):                                     # Loop thru k-mer range
        start_time = time.time()
        db_data.de_bruijn_graph(start=l_i, end=l_f, k=i)            # Create de Bruijn graph
        end_time = time.time()

        if(logging):
            file = db_data.make_docs(True, True, True, str(i))      # Create all files and get directed file
        else:
            file = db_data.directed_graph()                         # Gets dictionary of directed graph
    
        runtime.append([i, num_reads, end_time - start_time])

        if(res_align and logging):                                  # Perform alignment and create corresponding files
            file = eulerian_cycle_str(file, str(i), True, True)     # Create Eulerian cycle and path
            align_contig(file, i)
        elif(res_align and logging == False):                       # Perform alignment without corresponding files
            file = eulerian_cycle_str(file, str(i), False, False)   # Create Eulerian cycle and path
            align_contig(file, i)

    if(res_time):                                                   # Write program runtime
        utils.make_csv(runtime)
# ==============================================================================================================
def main():
    utils.rmove()
    file = "./input/sars_spike_protein_reads.fastq"
    loop_kmer(file, k_i=4, k_f=4, l_i=0, l_f=1, res_align=True, res_time=False, logging=False)
    

if __name__ == "__main__":
    main()