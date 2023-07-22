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

import matplotlib.pyplot as plt
import networkx as nx               # Plotting directed graph
import math
import copy
import random
import time

import utils
import alignment as al

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

        for s in seq:   # Iterate thru sequences with progress bar
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

        for k1 in kmers:
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
            for edge in self.edges:
                x1, x2 = edge
                f.write(x1 + '->' + x2 + '\n')
        f.close()

    # ----------------------------------------------------------------------------------------------------------
    # Writes edge data to a text file
    def create_directed_graph(self, file='output/temp/spike_protein_directed_graph.txt'):
        print("Creating Directed Graph File: ", file)

        with open(file, "w") as f:
            added_nodes = set()                                     # Set of (nodes, destination) pairs that have already been iterated through

            for edge in self.edges:
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
    # Create directed graph dictionary without writing to a file
    def directed_graph(self):
        dir_graph = {}
        added_nodes = set()                                     # Set of (nodes, destination) pairs that have already been iterated through

        for edge in self.edges:
            node,dest = edge
            if edge not in added_nodes:                         # Check if (node, destination) pair has already been writted
                dir_graph[node] = [dest]                        # Adding node and destination node to dictionary

                for edge2 in self.edges:                        # Iterate through all edges to see if there are any other edges coming out of node
                    node2, dest2 = edge2

                    # If the nodes are the same AND (node, destination) pair is not the same as above AND this node has not already been created
                    if node == node2 and edge != edge2 and edge2 not in added_nodes:
                        dir_graph[node] += [dest2]              # Adding additional destination for exising node
                        added_nodes.add(edge2)                  # add (node, destination) pair to already added set

        return dir_graph

    # ----------------------------------------------------------------------------------------------------------
    # Creates graph of connected nodes with corresponding kmers using matplotlib
    def matplot_graph(self, show_lab=False, show_fig=True, save_fig=False, file='./output/graph/deBruijn.png'):
        print("Creating Garph Image: ", file)
        
        plt.clf()
        fig = nx.DiGraph()                                              # Initialize weighted graph
        fig.add_edges_from(self.edges)                                  # Add edges to weighted graph
        k = 0.5/math.sqrt(fig.order())                                  # Used for spacing in spring layout
        #pos = nx.spring_layout(fig, k=k)
        pos = nx.shell_layout(fig, scale=2)                             # Set layout to shell (circular)

        options = {
            "node_color": "#A0CBE2",
            "node_size": 20,
            "edge_color": "#7d0901",
            "with_labels": show_lab,                                    # Show kmer labels
            "font_size": 5,
            "font_color": "#0a0a0a"
        }

        nx.draw(fig, pos, **options)                                    # Drawing directed graph
        plt.axis("off")                                                 # Do not show any axis
        if(show_fig): plt.show()                                        # Toggel show
        if(save_fig): plt.savefig(file, transparent=False ,dpi=500)     # Saving file and setting size

    # ----------------------------------------------------------------------------------------------------------
    # Master function for creating documents
    # edge_graph: creates graph image
    # edge_file: creates edge file used for directed graph
    # dir_graph: creates directed edge file
    # i: label for k-mer used (used for kmer loop)
    def make_docs(self, edge_graph=False, edge_file=False, dir_graph=False, i='1'):

        if(edge_graph): self.matplot_graph(True, False,True, './output/graph/deBruijn_' + i + '.png')
        if(edge_file): self.create_edges_file('output/temp/edges_' + i + '.txt')

        if(dir_graph): 
            file = 'output/temp/spike_protein_directed_graph_' + i + '.txt'
            self.create_directed_graph(file)
            return file

# ==============================================================================================================
# Eulerian Cycle
# ==============================================================================================================
def EulerianCycle(strings, format=True):
  ##Ignore this formatting block, it's only for a desired input
    if format:
        graph = [i.split(' -> ') for i in strings]
        graph = dict(graph)
        for (key, val) in graph.items():
            val = val.split(',')
            graph[key] = val

        copy_graph = copy.deepcopy(dict(graph))
    else:
        copy_graph = copy.deepcopy(strings)

    #We need to calculate the length of our graph so that we know when to stop the algorithm. This means we
    #need to calculate the number of edges in our graph.
    l = 1
    values = []
    for val in copy_graph.values():
        for i in val:
            l += 1
            values.append(i)

    #validate is set as default by true, you don't need this
    validate = True
    degrees = {key: [] for key in copy_graph.keys()}

    #Here, we're calculating the degree of the graph to determine whether each node has even degrees.
    for key in degrees.keys():
        degrees[key].append(values.count(key))
        degrees[key].append(len(copy_graph[key]))
        degrees[key].append(degrees[key][1] - degrees[key][0])

    #The final sequence is what will be our answer
    final_sequence = []

    #The cycle list will be the current cycle we are on
    cycle = []

    #Initialize a Eulerian cycle with a random node
    cn = random.choice([key for key in copy_graph.keys()])

    while len(final_sequence) != l:

        #The way we're going to pass through our graph is by removing nodes as we pass over all of their respective edges.
        #We keep on adding to the cycle as a result and the node we removed becomes our next node.

        if copy_graph[cn] != []:
            cycle.append(cn) #Add the start node to our cycle
            next_possibles = copy_graph[cn]
            new_cn = random.choice(next_possibles) #Using our dictionary, we find what nodes are connected to our current
                                                   #starting node and pick a random node to walk through next.
            #Remove the edge we've already walked through
            copy_graph[cn].remove(new_cn)
            cn = new_cn #our new node becomes our next starting node

        #A dead end of the graph is found if the node we are on has no more edges to walk through. Thus our cycle is complete.
        elif copy_graph[cn] == []:
            #We add this end node at the end of our final sequence.
            final_sequence.insert(0, cn)
            if len(cycle) == 0: #We're checking to see if we've passed through all possible cycles here
                # final_sequence.append(final_sequence[0])
                break

            #If there are still more cycles, we backtrack to the previous node and continue our walk from there
            else:
                cn = cycle[-1] #Our new starting node is the previous node -> this new node is later checked to make sure
                               #that it has unvisited edges, otherwise this is also definitely part of our final sequence
                               #and then we backtrack to the next previous node.
                #To prevent reusing this node, we remove it from the current cycle
                cycle.pop()

    return final_sequence

# ==============================================================================================================
# Eulerian Path
# ==============================================================================================================
def EulerianPath(strings, format=True):
    #Similar formatting for turning txt file to DeBruijn graph in python dict form
    if format:
        graph = [i.split(' -> ') for i in strings]
        graph = dict(graph)
        for (key, val) in graph.items():
            val = val.split(',')
            graph[key] = val
        copy_graph = copy.deepcopy(dict(graph))
    else:
        copy_graph = copy.deepcopy(strings)

    #Calculate the length the same way we did for the Eulerian cycle
    l = 1
    values = []
    for val in copy_graph.values():
        for i in val:
            l += 1
            values.append(i)

    validate = True

    #For a Eulerian path to be true, the end node must have a degree of 1, and it most cases the end node does
    #not connect to any other node with an edge, hence the graph is nearly balanced.

    end = [val for val in values if val not in copy_graph.keys()] #check if end node is connected to any other nodes

    #In some cases, this is not always true so we can check whether there is an end node or not
    if end == []:
        pass

    #If the condition is true, then we need to format our final dictionary such that we add a key in our dictionary
    #where the end node is assigned an empty list. If you remember the cycle loop in the Eulerian cycle code, this
    #empty list here will help us identify whether we have reached the end of the cycle or not.

    else:
        copy_graph[end[0]] = []
    degrees = {key: [] for key in copy_graph.keys()}

    #Regardless, we can verify our start and end nodes by checking their in-out degrees
    for key in degrees.keys():
        degrees[key].append(values.count(key))
        degrees[key].append(len(copy_graph[key]))
        degrees[key].append(degrees[key][1] - degrees[key][0])

    final_sequence = []
    cycle = []

    #cn or our current start node is determined if its degree is odd (or 1 in this case) TRY
    cn = -1
    try:
        cn = [key for key in degrees.keys() if degrees[key][2] == abs(1)][0]
    except Exception:
        pass
    finally:
        # test on other indexs
        cn = list(degrees.keys())[2:3]
        cn = str(cn[0])

    #We implement the exact same loop here as in the Eulerian cycle only with a difference where
    #we already know the starting node
    while len(final_sequence) != l:
        if copy_graph[cn] != []:
            cycle.append(cn)
            next_possibles = copy_graph[cn]
            new_cn = random.choice(next_possibles)
            copy_graph[cn].remove(new_cn)
            cn = new_cn
        elif copy_graph[cn] == []:
            final_sequence.insert(0, cn)

            #The length of our current cycle will be 0 or cycle=[] if there are no more nodes to backtrack to, therefore
            #the current node that we are on right now is the start node since it has no more outgoing edges nor incoming
            #ones that can allow us to backtrack from.

            if len(cycle) == 0:
                break
            else:
                cn = cycle[-1]
                cycle.pop()

    return final_sequence

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
# Finds eulerian cycle, writes the data to a text file
def eulerian_cycle_str(data, kmer, enter_file=False, logging=False):

    cycle = None
    if(enter_file):                                                 # Use file that was entered
        text = []
        print("Reading: " + data)
        with open(data, 'r') as f:                                  # Open directed graph file
            for line in f:
                text.append(line.strip())                           # Add the edges and strip the newline characters

        cycle = EulerianCycle(text)                                 # Find eulerian cycle
    else:
        cycle = EulerianCycle(data, format=False)


    if(logging):
        cycle_file = 'output/eulerian/eulerianCycle_' + kmer + '.txt'   # File that will be written to

        with open(cycle_file, 'w') as f:                                # Write eulerian cycle to file
            contig = cycle[0]
            for c in cycle[1:]:                                         # Assembling a contig
                contig += c[len(c)-1]                                   # Get last character and add it to contig
            f.write(''.join(contig))                                    # Write assembled contig to file

        return cycle_file                                               # Return file with assembled sequence

    else:
        contig = cycle[0]
        for c in cycle[1:]:                                             # Assembling a contig
            contig += c[len(c)-1]                                       # Get last character and add it to contig

        return ''.join(contig)                                          # Return assembled contig

# ==============================================================================================================
# Finds eulerian path, writes the data to a text file
def eulerian_path_str(data, kmer, enter_file=False, logging=False):

    cycle = None
    if(enter_file):                                                 # Use file that was entered
        text = []
        print("Reading: " + data)
        with open(data, 'r') as f:                                  # Open directed graph file
            for line in f:
                text.append(line.strip())                           # Add the edges and strip the newline characters

        cycle = EulerianCycle(text)                                 # Find eulerian cycle
    else:
        cycle = EulerianCycle(data, format=False)


    if(logging):
        cycle_file = 'output/eulerian/eulerianPath_' + kmer + '.txt'    # File that will be written to 

        with open(cycle_file, 'w') as f:                                # Write eulerian cycle to file
            contig = cycle[0]
            for c in cycle[1:]:                                         # Assembling a contig
                contig += c[len(c)-1]                                   # Get last character and add it to contig
            f.write(''.join(contig))                                    # Write assembled contig to file

        return cycle_file                                               # Return file with assembled sequence

    else:
        contig = cycle[0]
        for c in cycle[1:]:                                             # Assembling a contig
            contig += c[len(c)-1]                                       # Get last character and add it to contig

        return ''.join(contig)                                          # Return assembled contig

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
def main():
    utils.rmove()
    file = "./input/sars_spike_protein_reads.fastq"
    loop_kmer(file, k_i=4, k_f=4, l_i=0, l_f=1, res_align=True, res_time=False, logging=False)
    

if __name__ == "__main__":
    main()