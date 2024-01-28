import numpy as np
import networkx as nx
import alignment as al
import matplotlib.pyplot as plt
from typing import List

import matplotlib, math, copy, random, time, utils
matplotlib.use('agg')

LOGGER = utils.LOGGER

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
    def __init__(self, sequences:List[str]=[], header:List[str]=[], k:int=3, cut:int=1, cycle:bool=False):
        '''
        Initializes the de Bruijn graph
        
        Parameter(s):
            sequences (List[str]): a series of sequences being analyzed
            header (List[str]): a series of heder info for the corresponding sequence
            k (int, default=3): size of the k-mer, length of the substrings created from the main sequence
            cut (int, default=1): size of the prefix and suffix of the k-mer

        Output(s):
            None
        '''
        
        if cut > k:
            LOGGER.error(f"The prefix/suffix size cannot exceed the size of the k-mer: Cut={cut}, K={k}")
            raise utils.InvalidInput(f"The prefix/suffix size cannot exceed the size of the k-mer: Cut={cut}, K={k}")
        
        self.header = header
        self.sequences = sequences
        self.k = k
        self.cut = cut              # Number of chars removed to form a prefix/suffix

        self.kmers = {}             # Stores all the possible k-mer combinations from the input sequences
        self.edges = set()          # Stores all the possible edges connecting the k-mers
        self.dir_graph = {}         # Stores all the directed edges from the k-mers
        self.assem_contigs = []     # Stores a list of assembled contigs
        self.final_sequence = ''    # The final assembled sequence

        self.de_bruijn_graph(cycle=cycle)

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

        # End exceeds sequence list, change it to list length
        if(end > len(self.sequences) or end == -1):
            LOGGER.warning(f"The ending index exceeds the length of the sequence: {len(self.sequences)}, end: {end}")
            end = len(self.sequences)

        # Cannot exceed bounds of the list self.sequences
        if(start < 0 or start > end or start == end):
            LOGGER.error(f"Invalid Input: {start} is not less than {end} or greater than 0!")
            raise utils.InvalidInput(f"Invalid Input: {start} is not less than {end} or greater than 0!")

        # Get kmers from this list of sequences
        self.get_kmers(self.sequences[start:end], cycle)
        self.get_edges()
        self.get_directed_graph()
        self.get_eulerian_cycle()

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
        with k-1 mers as keys and lists of corresponding prefixs as values. Example, kmer = ACTG, prefix = ACT, 
        and suffix = CTG. The edge is then added if the prefix is in the suffix dictionary.
        
        Parameter(s): None
        
        Output(s): None
        '''

        kmer_suffixes = {}                          # Dictionary to store kmers with the same suffix

        # Group k-mers based on their suffix
        for kmer in self.kmers:
            prefix = kmer[:-self.cut]              # excluding the last character
            suffix = kmer[self.cut:]               # excluding the first character

            if prefix in kmer_suffixes:
                kmer_suffixes[prefix].append(suffix)
            else:
                kmer_suffixes[prefix] = [suffix]

        #print(f"Kmer Suffixes:\n{kmer_suffixes}")

        # Loop thru all the k-mers based on their prefix
        for kmer in self.kmers:
            prefix = kmer[:-self.cut]

            # If the prefix is in the dictionary, add the edge
            if prefix in kmer_suffixes:
                for suffix in kmer_suffixes[prefix]:
                    self.edges.add((prefix, suffix))

    # ----------------------------------------------------------------------------------------------------------
    def get_directed_graph(self):
        '''
        Creates a directed graph dictionary containing a node as a key and a list of destinations as the value
        
        Parameter(s): None
        
        Output(s): None
        '''

        for edge in self.edges:
            node, dest = edge   # set(node, destination)

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
        
        directed_graph = copy.deepcopy(self.dir_graph)
        num_edges = sum(len(val_list) for val_list in directed_graph.values()) + 1

        # Degree key: [incoming nodes, outgoing nodes]
        degrees = {key: 
                   [sum(val == key for val_list in directed_graph.values() for val in val_list),
                    len(directed_graph[key])] for key in directed_graph.keys()}

        # Create a list of nodes with no incoming edges
        starting_nodes = [key for key in degrees if degrees[key][0] == 0]

        # Iterate thru all the elements with no incoming nodes
        for current_node in starting_nodes:
            cycle = []
            contigs = []

            # Run until all the edges have been accounted for
            while len(contigs) != num_edges:

                # In-coming node has no connections
                if current_node not in directed_graph:
                    current_node = cycle[-1]
                    cycle.pop()
                    continue

                # Current node has out going edges and is not empty
                if directed_graph[current_node]:
                    cycle.append(current_node)
                    next_possibles = directed_graph[current_node]

                    # Randomly select a new node to iterate from
                    new_cn = random.choice(next_possibles)
                    directed_graph[current_node].remove(new_cn)
                    current_node = new_cn

                # Current node has no out going edges and is empty
                elif directed_graph[current_node] == []:
                    # Append the kmer to the contig list
                    contigs.insert(0, current_node)

                    if not cycle:
                        break
                    else:
                        current_node = cycle[-1]
                        cycle.pop()
            
            # Add the unassembled contigs to the list
            self.assemble_contigs(contigs)

    # ----------------------------------------------------------------------------------------------------------
    def assemble_contigs(self, contigs):
        '''
        Joins the unassembled contigs together into a main contig string and adds it to a contig list

        Parameter(s): None

        Output(s): None
        '''
        
        assembled_contig = contigs[0]                   # Start with first unassembled contig in the list
        num_chars = len(assembled_contig) - self.cut

        for contig in contigs[1:]:                      # Loop thru the remaining contigs
            assembled_contig += contig[num_chars:]      # Append the last few letter(s) to the final sequence

        self.assem_contigs.append(assembled_contig)

    # ----------------------------------------------------------------------------------------------------------
    def create_edges_file(self, filename:str='./temp/edges.txt'):
        '''
        Creates an edges file for testing
        
        Parameter(s):
            filename (str, default=./temp/edges.txt): file that stores the edge data

        Output(s):
            A path to the saved file containing the edge data
        '''

        try:
            LOGGER.info(f"Writing de Bruijn edges to: {filename}")

            with open(filename, 'w') as f:
                for edge in self.edges:
                    x1, x2 = edge
                    f.write(f"{x1}->{x2}\n")

            return filename
        
        except FileNotFoundError as e:
            LOGGER.error(f"File not found error when writing de Bruijn edges to {filename}: {str(e)}")
        except PermissionError as e:
            LOGGER.error(f"Permission error when writing de Bruijn edges to {filename}: {str(e)}")
        except Exception as e:
            LOGGER.error(f"An error occurred when writing de Bruijn edges to {filename}: {str(e)}")

    # ----------------------------------------------------------------------------------------------------------
    def create_directed_graph_file(self, filename='./temp/directed_graph.txt'):
        '''
        Creates a directed graph dictionary containing a node as a key and a list of destinations as the value.
        Writes the results to a text file. Used for testing
        
        Parameter(s): None
        
        Output(s):
            A path to the saved file containing nodes and their corresponding list of destinations.
        '''

        try:
            LOGGER.info(f"Writing de Bruijn directed edges to: {filename}")

            # Write the directed graph data to the file
            with open(filename, "w") as f:
                for node, destinations in self.dir_graph.items():
                    f.write(node + ' -> ' + ','.join(destinations) + '\n')

            return filename
        
        except FileNotFoundError as e:
            LOGGER.error(f"File not found error when writing de Bruijn directed edges to {filename}: {str(e)}")
        except PermissionError as e:
            LOGGER.error(f"Permission error when writing de Bruijn directed edges to {filename}: {str(e)}")
        except Exception as e:
            LOGGER.error(f"An error occured when writing de Bruijn directed edges to {filename}: {str(e)}")

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
        
        try:
            LOGGER.info(f"Saving de Bruijn garph figure to: {filename}")

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
            plt.savefig(filename, transparent=False ,dpi=500)  # Saving file and setting size

            return filename
        
        except FileNotFoundError as e:
            LOGGER.error(f"File not found error when saving de Bruijn graph figure to {filename}: {str(e)}")
        except PermissionError as e:
            LOGGER.error(f"Permission error when saving de Bruijn graph figure to {filename}: {str(e)}")
        except Exception as e:
            LOGGER.error(f"An error occurred when saving de Bruijn graph figure to {filename}: {str(e)}")

    # ----------------------------------------------------------------------------------------------------------
    def make_docs(self, edge_file:str=None, dir_graph_file:str=None, plot_file:str=None):
        '''
        Master function for creating multiple documents such as an edge file, directed edge file, 
        and a plot of the constructed graph.

        Parameter(s):
            edge_file (str, default=None): writes the edges of the graph to a file if not None
            dir_graph_file (str, default=None): writes the directed edges to a file if not None
            plot_file (str, default=None): saves the plotted De Bruijn graph figure to a file if not None
        
        Output(s):
            A list of paths to the saved files and their corresponding data (deBruijn plot, edges file, 
            and/or directed graph file).
        '''
        files = []

        # Write the graph edges to a file
        if(edge_file): 
            file = self.create_edges_file(filename=edge_file)
            files.append(file)

        # Write the directed graph edges to a file
        if(dir_graph_file): 
            file = self.create_directed_graph_file(filename=dir_graph_file)
            files.append(file)

        # Save the plotted graph figure
        if(plot_file): 
            file = self.plot_graph(show_label=True, filename=plot_file)
            files.append(file)

        return files

    # ----------------------------------------------------------------------------------------------------------
    def align_sequence(self, ref_data, filename:str="assembled_alignment.pdf"):
        '''
        Aligns assembled contig with reference genome, creates text file of alignment, and plots comparison. Used 
        for testing.

        Parameter(s):
            ref_data (List[List[str]]): a list composed of a list of sequences and a list of header info
            filename (str, defualt=assembled_alignment.pdf): name of the file where the alignment results are to be stored

        Output(s):
            A path to the file where the alignment comparison plot is saved
        '''

        a = al.Alignment(               # Initialize alignment
            ref=ref_data[0][0], 
            seq=self.final_sequence, 
            gap_pen=-2,
            match_point=1,
            match_pen=-1, 
            ignore=True
        )

        file = a.plot_compare(filename=filename)

        return file
    
    # ----------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Constructs a info string that emables the user to print the class attributes.

        Parameter(s): None

        Output(s):
            A string with all the classes atributes (k-mers, edges, directed graph, contigs, final sequence) and their counts
        '''
        return (
            f"K-mers ({len(self.kmers)}):\n{self.kmers}\n\n"
            f"Edges ({len(self.edges)}):\n{self.edges}\n\n"
            f"Directed Graph ({len(self.dir_graph)}):\n{self.dir_graph}\n\n"
            f"Assembled Contigs ({len(self.assem_contigs)}):\n{self.assem_contigs}\n\n"
            f"Final Sequence:\n{self.final_sequence}"
            )

# ==============================================================================================================
def loop_kmer(data, k_i:int=3, k_f:int=3, l_i:int=0, l_f:int=2, align:str=None, record_runtime:bool=False, get_files:bool=False):
    '''
    Test Function that runs through multiple k-mers. Loops thru a k-mer range, builds a set of kmers and edges
    NOTE: Used for testing!

    Parameter(s):
        data (List[List[str]]): a list composed of a list of sequences and a list of header info
        k_i (int, default=3): starting k-mer size
        k_f (int, default=3): ending k-mer size
        l_i (int, default=0): starting index for sequence reads
        l_f (int, default=2): ending index for sequence reads
        align (str, default=None): filename containing the reference sequence, get alignment if not None
        record_runtime (bool, default=False): record run-time data if True
        get_files (bool, default=False): create data files if True (k-mers, edges, etc)

    Ouput(s):
    '''
    runtime = []
    response = {'files':[]}
    num_reads = l_f - l_i                           # Number of reads being assembled

    if align is not None:
        ref_sequence = utils.get_data(align)

    for i in range(k_i, k_f):                       # Loop thru k-mer range
        start_time = time.time()
        graph = De_bruijn(data[0], data[1], i)      # Create de Bruijn graph
        end_time = time.time()
        runtime.append([i, num_reads, end_time - start_time])

        if get_files:                               # Get file information
            files = graph.make_docs(
                edge_file=True, 
                dir_graph=True, 
                plot_type='.png', 
                file_type='txt'
            )
            
            response["files"].append(files)

        if align is not None:                       # Get plot alignment file
            file = graph.align_sequence(
                ref_data=ref_sequence[0][0],
                k=i,
                file_type='.png'
            )

            response["files"].append(file)
    
    if(record_runtime):                             # Write program runtime
        file = utils.make_csv(runtime)
        response["files"].append(file)

    return response

# ==============================================================================================================
def main():
    # Testing File
    #data = utils.get_data(filename="./input/assembly_test.fastq")
    #graph = De_bruijn(sequences=data[0], header=['Testing'], k=11)
    #print(graph.final_sequence)

    word1 = 'pneumonoultramicroscopicsilicovolcanoconiosis'
    '''
    pneum
     neumo
      eumon
       umono
        ...
         iosis
    '''
    word1_fragments = [
        'pneum', 'neumo', 'eumon', 'umono', 'monou', 
        'onoul', 'noult', 'oultr', 'ultra', 'ltram', 
        'trami', 'ramic', 'amicr', 'micro', 'icros', 
        'crosc', 'rosco', 'oscop', 'scopi', 'copic', 
        'opics', 'picsi', 'icsil', 'csili', 'silic', 
        'ilico', 'licov', 'icovo', 'covol', 'ovolc', 
        'volca', 'olcan', 'lcano', 'canoc', 'anoco', 
        'nocon', 'oconi', 'conio', 'onios', 'niosi', 
        'iosis'
    ]

    word2 = 'hello world'

    word_graph = De_bruijn(sequences=[word1], header=['Testing Word'], k=5, cycle=False)
    print(word_graph)

if __name__ == "__main__":
    main()