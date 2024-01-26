# Phylogeny Inference

'''
Generate the distances between sequences using neighbor joining and build a phylogenic tree.
Gap in both sequences are treated as a similarity.

Input: fna data

Output:
1) genetic-distances.txt: A tab-delimited table of pairwise distances between all sequences.
2) edges.txt: A tab delimited. Each row describes an edge in the tree.
3) tree.txt: NEWICK format with all edge distances and with only tips named. For example, (A:0.1,B:0.2,(C:0.3,D:0.4):0.5)

'''
import sys, os, utils
import numpy as np
from Bio import Phylo
from typing import List
from io import StringIO
np.set_printoptions(threshold=sys.maxsize, precision=10, linewidth=np.inf)

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')
LOGGER = utils.LOGGER

PATH = os.path.dirname(os.path.abspath(__file__))

# ==============================================================================================================
# Building The Phylogeny Tree
# ==============================================================================================================
class Phylogeny_Tree:
    '''
    Stores an instance of a phylogeny tree.

    Parameter(s):
        data (dict): dictionary containing joined nodes and distances
        header (List[str]): list of header info for labels

    Output(s): None
    '''
    def __init__(self, data:dict, header:List[str]):
        self.edges = data
        self.labels = {}
        self.tree = None

        # Create dictionary for edge labels (exclude root)
        for i in range(len(header) + 1):
            self.labels[header[i-1]] = str(i)

        self.newick_tree()

    # ----------------------------------------------------------------------------------------------------------
    def newick_tree(self):
        '''
        Builds the Newick string used for creating the phylogeny tree figure.
        
        Parameter(s):
            None, uses the egdes stored in class data attribute.
            
        Output(s):
            None, builds class attribute tree using a Newick formatted string.
        '''

        def build_newick(parent):
            # Node is a taxa species and not an identifier
            if parent not in self.edges:
                return parent
            # Node is an identifier
            else:
                children = []
                for child_id, distance in self.edges[parent].items():
                    child_newick = build_newick(child_id)
                    children.append(f"{child_newick}:{distance}")

                return f"({','.join(children)})"

        self.tree = build_newick(self.edges["root"])

    # ----------------------------------------------------------------------------------------------------------
    def write_edges(self, filename:str="./temp/edges.txt"):
        '''
        Writes the edges of the phylogeny tree to an output file.

        Parameter(s):
            filename (str): name of the file that the edges are to be written to

        Output(s):
            A file containing the edges of the phylogeny tree.
        '''

        try:
            LOGGER.info(f"Writing phylogeny edges to: {filename}")

            def preorder_traversal(parent):
                result = []

                if parent in self.edges:
                    # Iterate thru both children (two or less)
                    for child, distance in self.edges[parent].items():
                        if child in self.labels:
                            child = self.labels[child]

                        result.append((f"{parent}\t{child}\t{distance}\n"))
                        result.extend(preorder_traversal(child))

                return result

            # Get preorder list of the edges so they can be printed
            data = preorder_traversal(self.edges["root"])

            with open(filename, 'w', newline='') as file:
                for line in data:
                    file.write(line)

            return filename
        
        except FileNotFoundError as e:
            LOGGER.error(f"File not found error when writing phylogeny edges to {filename}: {str(e)}")
        except PermissionError as e:
            LOGGER.error(f"Permission error when writing phylogeny edges to {filename}: {str(e)}")
        except Exception as e:
            LOGGER.error(f"An unexpected error occurred when writing phylogeny edges to {filename}: {str(e)}")

    # ----------------------------------------------------------------------------------------------------------
    def plot_ptree(self, filename:str="./temp/tree.pdf"):
        '''
        Creates a figure of a phylogenetic tree and saves it to a file.

        Parameter(s):
            filename (str): name of the file that the phylogeny tree is to be saved to

        Output(s):
            A file containing a figure of the phylogeny tree.
        '''

        try:
            LOGGER.info(f"Saving plotted phylogeny tree to: {filename}")

            tree = Phylo.read(StringIO(self.tree), "newick")

            plt.figure(figsize=(10,8), dpi=100)
            Phylo.draw(tree, do_show=False)
            plt.savefig(filename)

            plt.close()

            return filename
        
        except FileNotFoundError as e:
            LOGGER.error(f"File not found error when saving phylogeny tree figure to {filename}: {str(e)}")
        except PermissionError as e:
            LOGGER.error(f"Permission error when saving phylogeny tree figure to {filename}: {str(e)}")
        except Exception as e:
            LOGGER.error(f"An unexpected error occurred when saving phylogeny tree figure to {filename}: {str(e)}")

    # ----------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Overloads the string variable for printing the tree data
        
        Parameter(s): None
        
        Output(s):
            A formatted string comprised of the tree data
        '''

        return f"Data:\n{self.edges}\nTree:\n{self.tree}"

# ==============================================================================================================
# Building The Phylogeny Data
# ==============================================================================================================
class Phylogeny(Phylogeny_Tree):
    '''
    Stores an instance of the taxa sequences phylogeny and calculates the evolutionary distances between them.

    Parameter(s):
        sequences (List[str]): list of taxa sequences being evaluated
        header (List[str]): list of header info for labeling taxa

    Output(s): None
    '''
    def __init__(self, sequences:List[str], header:List[str]):

        self.sequences = sequences
        self.header = header
        self.dmatrix = None             # Distance matrix
        self.edges = {}                 # Edges of phylogenic tree

        self.nmatrix = None             # Distance matrix used for joining neighbors
        self.nheader = None             # Header list used for joining neighbors

        self.build_distance_matrix()    # Build distance matrix
        self.neighbor_joining()         # Connect neighboring sequences

        super().__init__(self.edges, self.header)

    # ----------------------------------------------------------------------------------------------------------
    def calc_distance(self, seq1:str, seq2:str):
        '''
        Calculates the percent dissimilarity between two sequences by dividing the number of mismatches 
        by the total number of bases.

        Parameter(s):
            seq1 (str): first sequence being compared
            seq2 (str): second sequence being compared

        Output(s)
            A float representing the percentage of dissimilarity between the two sequences.
        '''
        
        match = 0
        mismatch = 0

        if(len(seq1) == 0 or len(seq2) == 0):           # Sequences cannot be null
            raise utils.InvalidInput("ERROR: NULL STRING!")
        if(len(seq1) != len(seq2)):                     # Lengths of comparing sequences are not equal
            raise utils.InvalidInput("ERROR: lengths are not equal!")

        total = len(seq1)                               # Total number of bases
        for (i,j) in zip(seq1,seq2):
            if(i == j):
                match += 1
            else:
                mismatch += 1

        return (mismatch/total)                         # Return % dissimilarity

    # ----------------------------------------------------------------------------------------------------------
    def build_distance_matrix(self):
        '''
        Builds the distance matrix of the using the input sequences. The input sequences must have equal length so 
        a nxn matrix can be created. Since the matrix is a reflection along the diagonal, only half of it is iterated 
        through and the values are reflected to their corresponding coordinates [i,j] = [j,i].

        Parameter(s): None

        Output(s): None
        '''

        size = len(self.sequences)                                      # Size of n x n distance matrix
        self.dmatrix = np.array([[0.0]*size for i in range(size)])      # Initialize distance matrix

        for i in range(size):                                           # Sequence 1
            for j in range(i, size):                                    # Sequence 2 (Compare)
                if(i == j): continue
                distance = self.calc_distance(self.sequences[i],self.sequences[j])        # Getting % dissimilarity between sequences
                #print(i , ", ", j, ", ", dis)
                self.dmatrix[i][j] = self.dmatrix[j][i] = distance      # Setting % dissimilarity in distance matrix

    # ----------------------------------------------------------------------------------------------------------
    def q_matrix(self):
        '''
        Calculates the Q matrix used for determining which neighbors to join. The input sequences must have equal length so 
        a nxn matrix can be created. Since the matrix is a reflection along the diagonal, only half of it is iterated 
        through and the values are reflected to their corresponding coordinates [i,j] = [j,i].

        Parameter(s): None

        Output(s):
            A Q matrix used for connecting neighboring taxa.
        '''
        # Invalid matrix demensions
        if(self.nmatrix.shape[0] != self.nmatrix.shape[1]):
            raise utils.InvalidInput(
                f"ERROR: NOT A N X N MATRIX!\nMatrix Shape: {self.nmatrix.shape[0]} X {self.nmatrix.shape[1]}"
                )

        size = self.nmatrix.shape[1]
        qmatrix = np.array([[0.0]*size for i in range(size)])   # Initialize Q matrix

        for i in range(0, size):                                # Iterate thru rows
            for j in range(i+1, size):                          # Iterate thru cols
                if(i == j): continue                            # Skip diagonal elements

                a_sum = 0.0
                for x in range(size):                           # Sum the values in a
                    a_sum += self.nmatrix[x][j]

                b_sum = 0.0
                for y in range(size):                           # Sum the values in b
                    b_sum += self.nmatrix[i][y]

                # Calculating Q values
                value = (self.nmatrix[i][j])*(size - 2) - a_sum - b_sum
                qmatrix[j][i], qmatrix[i][j] = value, value

        return qmatrix

    # ----------------------------------------------------------------------------------------------------------
    def get_distance(self, a:int, b:int):
        '''
        Calculates distances from each node in the neighbor distance matrix.

        Parameter(s):
            a (int): index location for first set of distances
            b (int): index location for second set of distances

        Output(s):
            The evolutionary distance between two taxa.
        '''

        if(a >= self.nmatrix.shape[0] or b >= self.nmatrix.shape[0]):
            raise utils.InvalidInput("ERROR: INPUT OUT OF RANGE OF MATRIX!")

        sum_a = sum(self.nmatrix[a][q] for q in range(self.nmatrix.shape[0]))   # Summing distances in a
        sum_b = sum(self.nmatrix[b][q] for q in range(self.nmatrix.shape[0]))   # Summing distances in b
        diff = sum_a - sum_b                                                    # Computing difference between a and b
        
        s = 2 * (self.nmatrix.shape[0] - 2)
        distance = (0.5)*(self.nmatrix[a][b]) + (diff / s)                      # Calculating distance from a to u

        return distance

    # ----------------------------------------------------------------------------------------------------------
    def join_neighbor(self, node:int):
        '''
        Joins neighboring sequences based on evolutionary distance. Updates a cached distance matrix and creates 
        a branch between the most related taxa. Used as a helper function for neighbor_joining.

        Parameter(s):
            node (int): a numberical label for marking taxa position

        Output(s): 
            None, manipulates the class variables nmatrix, nheader, and edges.
        '''

        if(self.nmatrix.shape[0] != self.nmatrix.shape[1]):
            raise utils.InvalidInput(
                f"ERROR: MATRIX MUST BE N X N! Matix shape is {self.nmatrix.shape[0]} X {self.nmatrix.shape[1]}."
            )

        qmatrix = self.q_matrix()               # Build Q matrix
        min = np.amin(qmatrix)                  # Min value
        loc = np.where(qmatrix == min)          # Find mins

        a, b = loc[0][0], loc[0][1]             # Index location of first min value
        d_au = self.get_distance(a, b)          # Calculating distance from a to u
        d_bu = self.nmatrix[a][b] - d_au        # Calculating distance from b to u

        #print(f"Min: {qmatrix[a][b]} d_ab: {self.nmatrix[a][b]} d_au: {d_au} d_bu: {d_bu}\n")

        u = []
        for i in range(self.nmatrix.shape[0]):
            if(i == a): continue
            # Calaculating distances of each taxa to u
            d_u = (0.5)*(self.nmatrix[a][i] + self.nmatrix[i][b] - self.nmatrix[a][b])
            u.append(d_u)

        self.nmatrix = np.delete(self.nmatrix, [a], axis=0)                 # Deleting row a
        self.nmatrix = np.delete(self.nmatrix, [a], axis=1)                 # Deleting column a

        self.nmatrix[b-1, 0:self.nmatrix.shape[0]] = u                      # Adding u distances to row (b moved 1)
        self.nmatrix[0:self.nmatrix.shape[0], b-1] = u                      # Adding u distances to column (b moved 1)

        self.edges[node] = {self.nheader[a]: d_au, self.nheader[b]: d_bu}   # Branch from a to u
        self.nheader[b] = node                                              # Add new node to list
        self.nheader.pop(a)                                                 # Remove extra label

    # ----------------------------------------------------------------------------------------------------------
    # Builds list of branche edges
    def neighbor_joining(self):
        '''
        Joins all the neighboring taxa together based on the calculated evolutionary distance. Stores the results 
        in self.edges used for creating the phylogeny tree.

        Parameter(s): None

        Output(s): 
            None, results are stored in self.edges for tree construction.
        '''

        if self.dmatrix is None:
            self.build_distance_matrix()
        
        self.nmatrix = self.dmatrix.copy()
        self.nheader = self.header.copy()
        node = len(self.header) + 1             # Joining node label

        while(self.nmatrix.shape[0] > 3):       # Iterate thru taxa until there are only two nodes left
            self.join_neighbor(node)            # Join a neighboring taxa
            node += 1                           # Updating new node label

        d_vw = self.get_distance(0, 1)          # Distance from v to w
        d_wd = self.nmatrix[0][1] - d_vw        # Distance from d to w
        d_we = self.nmatrix[0][2] - d_vw        # Distance from e to w

        # Add last three nodes
        self.edges[node] = {self.nheader[0]:d_vw, self.nheader[1]:d_wd, self.nheader[2]:d_we}
        self.edges['root'] = node               # Track root node

    # ----------------------------------------------------------------------------------------------------------
    def write_distance_matrix(self, filename:str="./temp/genetic-distances.txt"):
        '''
        Writes the distance matrix to a file.
        
        Parameter(s):
            filename (str, default=./temp/genetic-distances.txt): the file that the distance matrix will be written to

        Output(s):
            A path to the file containing the distance matrix data.            
        '''

        try:
            LOGGER.info(f"Writing phylogeny distance matrix to: {filename}")

            if not self.header:
                raise ValueError("Header is not defined.")

            if self.dmatrix is None:
                raise ValueError("Distance matrix is not defined.")

            with open(filename, 'w') as file:
                file.write('\t'.join(self.header) + '\n')               # Write column names from header

                for x in range(len(self.header)):                       # Iterate thru rows in the distance matrix
                    s = '\t'.join([str(i) for i in self.dmatrix[x]])    # Join the row data together
                    file.write(f"{self.header[x]}\t{s}\n")              # Write the header + data for each corresponding row

            return filename
        
        except FileNotFoundError as e:
            LOGGER.error(f"File not found error when writing phylogeny distance matrix to {filename}: {str(e)}")
        except PermissionError as e:
            LOGGER.error(f"Permission error when writing phylogeny distance matrix to {filename}: {str(e)}")
        except Exception as e:
            LOGGER.error(f"An unexpected error occurred when writing phylogeny distance matrix to {filename}: {str(e)}")

    # ----------------------------------------------------------------------------------------------------------
    def get_files(self, 
            matrix_file:str="./temp/genetic-distances.txt",
            edge_file:str="./temp/edges.txt",
            tree_file:str="./temp/tree.pdf"
        ):
        '''
        Get all the phylogeny files used for constructing the tree.

        Parameter(s):
            matrix_file (str, default=./temp/genetic-distances.txt): file name where the distance matrix data is stored
            edge_file (str, default=./temp/edges.txt): file name where the edge data is stored
            tree_file (str, default=./temp/tree.pdf): file name where the plotted tree figure is stored

        Output(s):
            A list of files containing the distance matrix, tree edges, and plotted phylogeny tree.
        '''
        files = []
    
        files.append(self.write_distance_matrix(filename=matrix_file))
        files.append(self.write_edges(filename=edge_file))
        files.append(self.plot_ptree(filename=tree_file))
    
        return files
    
    # ----------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Overloads the string operator to display class data.
        
        Parameter(s): None
        
        Output(s):
            A string composed of the header info, initial distance matrix, and tree edges (neighbors).
        '''
        return (
            f"Header:\n{self.header}\n\n"
            f"Distance Matrix:\n{self.dmatrix}\n\n"
            f"Edges:\n{self.edges}\n\n"
            f"Newick Tree:\n{self.tree}"
        )

# ================================================================================================================================================
if __name__ == "__main__":
    data = utils.get_data("./input/phylogeny_test2.fna")

    phyl = Phylogeny(data[0], data[1])
    phyl.get_files()