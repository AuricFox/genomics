# Phylogeny Inference

'''
Generate the distances between sequences using neighbor joining and build a phylogenic tree.
Gap in both sequences are treated as a similarity.

Input: fna data

Output:
1) genetic-distances.txt: A tab-delimited table of pairwise distances between all sequences.
2) edges.txt: A tab delimited. Each row describes an edge in the tree.
3) tree.txt: NEWICK format with all edge distances and with only tips named. For example, (A:0.1,B:0.2,(C:0.3,D:0.4):0.5)

NOTE: Node class is an adaption of the TreeNode class from scikit-bio.
https://github.com/biocore/scikit-bio/blob/0.2.3/skbio/tree/_tree.py#L53

'''
import sys, os, utils
import numpy as np
from typing import List
from collections import defaultdict
np.set_printoptions(threshold=sys.maxsize, precision=10, linewidth=np.inf)

# ==============================================================================================================
# Tree Nodes
# ==============================================================================================================
class Node:
    '''
    Stores an instance of a tree data structure.

    Parameter(s):
        name (str, default=None): name of the species or type.
        distance (float, default=None): evolutionary distance between species/nodes.
        parent (Node, default=None): parent node.
        children (List[Node], default=None): nodes of existing children.

    Output(s): None
    '''
    def __init__(self, name:str=None, distance:float=None, parent:'Node'=None, children:List['Node']=None):
        self.name = name
        self.distance = distance
        self.parent = parent

        self._tip_cache = {}
        self._non_tip_cache = {}
        self._registered_caches = set()
        self.children = []

        if children is not None:
            self.extend(children)

    # ----------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Adds children to the tree.
        
        Parameter(s): None
            
        Output(s):
            A string of the tree data comprised of names and distances
        '''
        return self.name
    
    # ----------------------------------------------------------------------------------------------------------
    def is_root(self):
        '''
        Returns True if the current node is the root.
        
        Parameter(s): None
        
        Output(s):
            True if the current node is the root, else False if the node is not the root.
        '''

        return self.parent is None
    
    # ----------------------------------------------------------------------------------------------------------
    def root(self):
        '''
        Gets the root of the tree.
        
        Parameter(s): None
        
        Output(s): None
        '''

        node = self
        while not node.is_root():
            node = node.parent

        return node
    
    # ----------------------------------------------------------------------------------------------------------
    def invalidate_caches(self, attr=True):
        
        if not self.is_root():
            self.root().invalidate_caches()
        else:
            self._tip_cache = {}
            self._non_tip_cache = {}

            if self._registered_caches and attr:
                for n in self.traverse():
                    for cache in self._registered_caches:
                        if hasattr(n, cache):
                            delattr(n, cache)

    # ----------------------------------------------------------------------------------------------------------
    def create_caches(self):
        '''
        Constructs a lookup cache for node names.
        
        Parameter(s): None
        
        Output(s): None
        '''

        if not self.is_root():
            self.root().create_caches()
        else:
            if self._tip_cache and self._non_tip_cache:
                return

            self.invalidate_caches(attr=False)

            tip_cache = {}
            non_tip_cache = defaultdict(list)

            for node in self.postorder():
                name = node.name

                if name is None:
                    continue

                if node.is_tip():
                    if name in tip_cache:
                        raise utils.InvalidInput(f"Tip with name {name} already exists!")

                    tip_cache[name] = node
                else:
                    non_tip_cache[name].append(node)

            self._tip_cache = tip_cache
            self._non_tip_cache = non_tip_cache

    # ----------------------------------------------------------------------------------------------------------
    def _adopt(self, node:'Node'):
        self.invalidate_caches()

        if node.parent is not None:
            node.parent.remove(node)
        node.parent = self

        return node
    # ----------------------------------------------------------------------------------------------------------
    def extend(self, nodes:List['Node']):
        '''
        Adds children to the tree.
        
        Parameter(s):
            nodes (List[Node]): a list of children to be added to the tree.
            
        Output(s): None
        '''

        self.children.extend([self._adopt(n) for n in nodes])
    
    # ----------------------------------------------------------------------------------------------------------
    def find(self, name:str):
        '''
        Searches for the first instance of a node with the same name.
        
        Parameter(s):
            name (str): name of the node being searched.
            
        Output(s): None
        '''
        root = self.root()

        # if what is being passed in looks like a node, just return it
        if isinstance(name, root.__class__):
            return name

        root.create_caches()
        node = root._tip_cache.get(name, None)

        if node is None:
            node = root._non_tip_cache.get(name, [None])[0]

        if node is None:
            raise utils.InvalidInput(f"Node {name} is not in self!")
        else:
            return node
    
    # ----------------------------------------------------------------------------------------------------------
    def append(self, node:'Node'):
        '''
        Appends a node to children.

        Parameter(s):
            node (Node): node being appended to list of children
        
        Output(s): None
        '''
        self.children.append(self._adopt(node))

    # ----------------------------------------------------------------------------------------------------------
    def pop(self, index:int=-1):
        '''
        Remove a Node from the tree.

        Parameter(s):
            idx (int): index location of the node being removed.
        
        Output(s):
            The node removed from the tree. 
        '''
        return self._remove_node(index)

    # ----------------------------------------------------------------------------------------------------------
    def _remove_node(self, index:int):
        '''
        Performs node removal.
        
        Parameter(s):
            idx (int): index location of the node being removed.
        
        Output(s):
            The node removed from the tree.        
        '''
        self.invalidate_caches()
        node = self.children.pop(index)
        node.parent = None

        return node

    # ----------------------------------------------------------------------------------------------------------
    def remove(self, node:'Node'):
        '''
        Removes a node from the tree.

        Parameter(s):
            node (Node): the current node being removed.

        Output(s):
            True if the node was removes, else returns false.
        '''
        for (i, curr_node) in enumerate(self.children):
            if curr_node is node:
                self._remove_node(i)
                return True
            
        return False
    
    # ----------------------------------------------------------------------------------------------------------
    def traverse(self, self_before=True, self_after=False, include_self=True):
        '''
        Iterates over descendants and returns them.
        
        Paramter(s):
            self_before (bool, default=True): includes each node before its descendants if True.
            self_after (bool, default=False): includes each node after its descendants if True.
            include_self (bool, default=False): include the initial node if True.
            
        Output(s):
            Yields successive Node objects.
        '''
        if self_before:
            if self_after:
                return self.pre_and_postorder(include_self=include_self)
            else:
                return self.preorder(include_self=include_self)
        else:
            if self_after:
                return self.postorder(include_self=include_self)
            else:
                return self.tips(include_self=include_self)

    # ----------------------------------------------------------------------------------------------------------
    def preorder(self, include_self=True):
        '''
        Performs preorder iteration over the tree.
        
        Parameter(s):
            include_self (bool, default=True): include the initial node if True.
            
        Output(s):
            Yields successive Node objects.
        '''
        
        stack = [self]
        while stack:
            curr = stack.pop()
            if include_self or (curr is not self):
                yield curr
            if curr.children:
                stack.extend(curr.children[::-1])

    # ----------------------------------------------------------------------------------------------------------
    def postorder(self, include_self=True):
        '''
        Performs postorder iteration over tree.
        
        Parameter(s):
            include_self (bool, default=True): include the initial node if True.

        Output(s):
            Yields successive Node objects.
        '''
        child_index_stack = [0]
        curr = self
        curr_children = self.children
        curr_children_len = len(curr_children)
        while 1:
            curr_index = child_index_stack[-1]
            # if there are children left, process them
            if curr_index < curr_children_len:
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_children_len = len(curr_children)
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                curr_children_len = len(curr_children)
                child_index_stack.pop()
                child_index_stack[-1] += 1

    # ----------------------------------------------------------------------------------------------------------
    def pre_and_postorder(self, include_self=True):
        '''
        Performs iteration over tree, visiting node before and after.
        
        Parameter(s):
            include_self (bool, default=True): include the initial node if True.

        Output(s):
            Yields successive Node objects.
        '''
        # handle simple case first
        if not self.children:
            if include_self:
                yield self
            raise StopIteration
        child_index_stack = [0]
        curr = self
        curr_children = self.children
        while 1:
            curr_index = child_index_stack[-1]
            if not curr_index:
                if include_self or (curr is not self):
                    yield curr
            # if there are children left, process them
            if curr_index < len(curr_children):
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                child_index_stack.pop()
                child_index_stack[-1] += 1

    # ----------------------------------------------------------------------------------------------------------
    def is_tip(self):
        '''
        Finds if the node is a tip (has no children) or not.

        Parameter(s): None

        Output(s):
            Returns True if the node is a tip and has no children, else False.
        '''
        return not self.children
    
    # ----------------------------------------------------------------------------------------------------------
    def tips(self, include_self=False):
        '''
        Iterates over tips descended from self.

        Parameter(s):
            include_self (bool, default=False): include the initial node if True

        Output(s):
            Yeilds successive Node objects.
        '''
        for node in self.postorder(include_self):
            if node.is_tip():
                yield node

    # ----------------------------------------------------------------------------------------------------------
    def non_tips(self, include_self=False):
        '''
        Iterates over nontips descended from self.

        Parameter(s):
            include_self (bool, default=False): include the initial node if True

        Output(s):
            Yeilds successive Node objects.
        '''
        for node in self.postorder(include_self):
            if not node.is_tip():
                yield node

# ==============================================================================================================
# Building The Phylogeny Tree
# ==============================================================================================================
class Ptree:
    '''
    Stores an instance of a phylogeny tree.

    Parameter(s):
        data (dict): dictionary containing joined nodes and distances
        header (List[str]): list of header info for labels

    Output(s): None
    '''
    def __init__(self, data:dict, header:List[str]):
        self.data = data
        self.tree = None                    # Tree data structure
        self.header = {}

        # Create dictionary for edge labels (exclude root)
        for i in range(len(header) + 1):
            self.header[header[i-1]] = str(i)

        self.init_tree()

    # ----------------------------------------------------------------------------------------------------------
    def init_tree(self):
        '''
        Initializes construction of the phylogeny tree and calls build_tree that recursively adds nodes to the tree.

        Parameter(s): None

        Output(s): None
        '''

        node = self.data['root']
        self.tree = Node(str(node))                   # Set root node

        keys = []                                   # List of keys for connected nodes
        for k in self.data[node]:                   # Get all the keys (should be 3)
            keys.append(k)

        if(len(keys) != 3):                         # Root can only have 3 children
            print("ERROR: INCORRECT NUMBER OF CHILDREN!")
            return

        child1 = Node(str(keys[0]), self.data[node][keys[0]], parent=self.tree)
        child2 = Node(str(keys[1]), self.data[node][keys[1]], parent=self.tree)
        child3 = Node(str(keys[2]), self.data[node][keys[2]], parent=self.tree)
        self.tree.extend([child1, child2, child3])  # Add Children to tree (root)

        self.build_tree(keys[0])                     # Build off child 1
        self.build_tree(keys[1])                     # Build off child 2
        self.build_tree(keys[2])                     # Build off child 3

    # ----------------------------------------------------------------------------------------------------------
    def build_tree(self, node:'Node'):
        '''
        Recursively builds the phylogeny tree by adding nodes to the branches.
        
        Parameter(s)
            node (Node): current node being added to the tree
            
        Output(s): None
        '''
        
        if node in self.data:           # Is node in the dictionary
            keys = []                   # List of keys for connected nodes
            for k in self.data[node]:   # Get all the keys (should be 2)
                keys.append(k)

            if(len(keys) != 2):         # Can only have two children
                print("ERROR: INCORRECT NUMBER OF CHILDREN!")
                return

            # Children of the input node
            child1 = Node(str(keys[0]), self.data[node][keys[0]])   # Creating child 1
            child2 = Node(str(keys[1]), self.data[node][keys[1]])   # Creating child 2
            self.tree.find(str(node)).extend([child1,child2])       # Find parent and add children

            self.build_tree(keys[0])                                # Build off child 1
            self.build_tree(keys[1])                                # Build off child 2

        return

    # ----------------------------------------------------------------------------------------------------------
    # Writes tree data to output file
    def write_tree(self, filename="./temp/tree.tre"):
        self.tree.write(filename)

    # ----------------------------------------------------------------------------------------------------------
    def write_edges(self, filename:str="./temp/edges.txt"):
        '''
        Writes the edges of the phylogeny tree to an output file.

        Parameter(s):
            filename (str): name of the file that the edges are to be written to

        Output(s):
            A file containing the edges of the phylogeny tree.
        '''

        with open(filename, 'w', newline='') as file:

            for x in self.tree.preorder():
                if(x.distance == None): continue            # Skip root node

                name = x.name
                if(name in self.header):
                    name = self.header[x.name]              # Convert to index

                line =  f"{x.parent.name}\t{name}\t{str(x.distance)}\n"
                file.write(line)

    # ----------------------------------------------------------------------------------------------------------
    def create_tree(self, filename:str="./temp/tree.pdf"):
        '''
        Runs R script to create a pdf of the phylogeny tree.

        Parameter(s):
            filename (str): name of the file that the phylogeny tree is to be saved to

        Output(s):
            A file containing a figure of the phylogeny tree.
        '''

        input = f"Rscript hw3-plot-edges.r ./output/edges.txt hw3-tip-labels.txt {filename}"
        os.system(input)

    # ----------------------------------------------------------------------------------------------------------
    def creat_newick(self, filename:str="./temp/tree-newick.pdf"):
        '''
        Runs R script to create a pdf of the newick version of a phylogeny tree.

        Parameter(s):
            filename (str): name of the file that the phylogeny tree is to be saved to

        Output(s):
            A file containing a figure of the newick tree.
        '''

        input = f"Rscript hw3-plot-newick.r ./output/tree.tre hw3-tip-labels.txt {filename}"
        os.system(input)

    # ----------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Overloads the string variable for printing the tree data
        
        Parameter(s): None
        
        Output(s):
            A formatted string comprised of the tree data
        '''

        return f"Data:\n{self.data}\nTree:\n{self.tree}"

# ==============================================================================================================
# Building The Phylogeny Data
# ==============================================================================================================
class Phylogeny:
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
        self.tree = None                # Nodes of phylogenic tree

        self.nmatrix = None             # Distance matrix used for joining neighbors
        self.nheader = None             # Header list used for joining neighbors

        self.build_distance_matrix()    # Build distance matrix
        self.neighbor_joining()         # Connect neighboring sequences

    # ----------------------------------------------------------------------------------------------------------
    def calc_distance(self, seq1, seq2):
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

        Output(s): None
        '''
        # Invalid matrix demensions
        if(self.nmatrix.shape[0] != self.nmatrix.shape[1]):
            raise utils.InvalidInput(
                f"ERROR: NOT A N X N MATRIX!\nMatrix Shape: {self.nmatrix.shape[0]} X {self.nmatrix.shape[1]}"
                )

        size = self.nmatrix.shape[1]
        qmatrix = np.array([[0.0]*size for i in range(size)])               # Initialize Q matrix

        for i in range(0, size):                                            # Iterate thru rows
            for j in range(i+1, size):                                      # Iterate thru cols
                if(i == j): continue                                        # Skip diagonal elements

                a_sum = 0.0
                for x in range(size):                                       # Sum the values in a
                    a_sum += self.nmatrix[x][j]

                b_sum = 0.0
                for y in range(size):                                       # Sum the values in b
                    b_sum += self.nmatrix[i][y]

                qmatrix[j][i] = (self.nmatrix[i][j])*(size - 2) - a_sum - b_sum # Calculating Q values
                qmatrix[i][j] = (self.nmatrix[i][j])*(size - 2) - a_sum - b_sum # Calculating Q values

        return qmatrix

    # ----------------------------------------------------------------------------------------------------------
    def get_distance(self, a:int, b:int):
        '''
        Calculates distances from each node in the neighbor distance matrix.

        Parameter(s):
            a (int): index location for first set of distances
            b (int): index location for second set of distances

        Output(s):
            The
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
    # Joins neighboring sequences
    # Returns a tuple containing (new dmatrix, new list of header info, edge A, edge B)
    def join_neighbor(self, node):

        if(self.nmatrix.shape[0] != self.nmatrix.shape[1]):
            print("ERROR: MATRIX MUST BE N X N!")
            return

        qm = self.q_matrix(self.nmatrix)                # Build Q matrix
        min = np.amin(qm)                               # Min value
        loc = np.where(qm == min)                       # Find mins

        #print('\n', loc)
        #print(m,'\n\n',qm,'\n')

        a, b = loc[0][0], loc[0][1]                     # Index location of first min value
        d_au = self.get_distance(a, b)                  # Calculating distance from a to u
        d_bu = self.nmatrix[a][b] - d_au                # Calculating distance from b to u
        #print("Min: ", qm[a][b], " d_ab: ", m[a][b], " d_au: ", d_au, " d_bu: ", d_bu, '\n')

        u = []
        for i in range(self.nmatrix.shape[0]):
            if(i == a): continue
            # Calaculating distances of each taxa to u
            d_u = (0.5)*(self.nmatrix[a][i] + self.nmatrix[i][b] - self.nmatrix[a][b])
            u.append(d_u)

        self.nmatrix = np.delete(self.nmatrix, [a], axis=0)                 # Deleting row a
        self.nmatrix = np.delete(self.nmatrix, [a], axis=1)                 # Deleting column a
        #print("Trimming: \n", m, "\n")

        self.nmatrix[b-1, 0:self.nmatrix.shape[0]] = u                      # Adding u distances to row (b moved 1)
        self.nmatrix[0:self.nmatrix.shape[0], b-1] = u                      # Adding u distances to column (b moved 1)
        #print("Adding new distances: ", u, "\n", m)

        self.edges[node] = {self.nheader[a]: d_au, self.nheader[b]: d_bu}   # Branch from a to u
        self.nheader[b] = node                                              # Add new node to list
        self.nheader.pop(a)                                                 # Remove extra label

    # ----------------------------------------------------------------------------------------------------------
    # Joins all the neighboring taxa together
    # Builds list of branche edges
    def neighbor_joining(self):

        if self.dmatrix is None:
            self.build_distance_matrix()

        
        self.nmatrix = self.dmatrix.copy()
        self.nheader = self.header.copy()
        node = len(self.header) + 1             # Joining node label

        while(self.nmatrix.shape[0] > 3):       # Iterate thru taxa until there are only two nodes left
            data = self.join_neighbor(node)     # Join a neighboring taxa
            self.nmatrix = data[0]              # Updating matrix
            self.nheader = data[1]              # Updating header
            node += 1                           # Updating new node label

            #print(f"Distance Matrix:\n{dmatrix}\nJoining: {data[2:4]}\n")

        #print(self.q_matrix(m))
        d_vw = self.get_distance(0, 1)          # Distance from v to w
        d_wd = self.nmatrix[0][1] - d_vw        # Distance from d to w
        d_we = self.nmatrix[0][2] - d_vw        # Distance from e to w

        # Add last three nodes
        self.edges[node] = {self.nheader[0]:d_vw, self.nheader[1]:d_wd, self.nheader[2]:d_we}
        self.edges['root'] = node               # Track root node

        #print(f"Distance v-w: {d_vw}\nDistance w-d: {d_wd}\nDistance w-e: {d_we})

    # ----------------------------------------------------------------------------------------------------------
    # Prints out values for debugging
    def debug(self):
        #print(self.header, '\n\n', self.dmatrix, '\n\n')
        print(self.edges)

# ==============================================================================================================
# Utility Functions
# ==============================================================================================================
# Retreives data from the fna file and returns a tuple containing a list of sequences and headers
def get_data(filename):
    file = open(filename, 'r')
    seq_data = []
    header_data = []

    for line in file:                                   # Read each line in file
        line = line.strip()                             # Strip newline characters

        if(line == ""):                                 # Empty line
            # print("Not text")
            continue
        elif(line[0] == ">"):                           # Header information
            header_data.append(line[1:])
        else:                                           # Genetic sequence
            seq_data.append(line)

    file.close()
    return (seq_data, header_data)

# ==============================================================================================================
# Write distance matrix to a text file
def write_dmatrix(data, filename="./output/genetic-distances.txt"):

    with open(filename, 'w', newline='') as file:
        file.write('\t'.join(data.header) + '\n')               # Write column names from header

        for x in range(len(data.header)):                       # Iterate thru rows in the distance matrix
            s = '\t'.join([str(i) for i in data.dmatrix[x]])    # Join the row data together
            file.write(data.header[x] + '\t' + s + '\n')        # Write the header + data for each corresponding row

# ==============================================================================================================
# Write tree matrix to a text file
def write_data(data, header, filename="./output/edges.txt"):
    tr = Ptree(data, header)        # Initializing tree
    tr.write_edges()                # Writing edges to file
    tr.write_tree()                 # Writing newick tree to file

# ==============================================================================================================
def main():
    return

# ================================================================================================================================================
if __name__ == "__main__":
    main()
