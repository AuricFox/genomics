# Samuel Kitzerow, kitze012
# Phylogeny Inference
# Creates the edge and tree data structures

from skbio import TreeNode as tn
import os

class Ptree:
    def __init__(self, data, header):
        self.data = data        # Dictionary containing joined nodes and distances
        self.tree = None        # Tree data structure
        self.header = {}

        for i in range(len(header) + 1):    # Create dictionary for edge labels
            self.header[header[i-1]] = str(i)

        self.tree_data()

    # ======================================================================
    # Builds the tree recusively
    def help_tree(self, node):
        if node in self.data:           # Is node in the dictionary
            keys = []                   # List of keys for connected nodes
            for k in self.data[node]:   # Get all the keys (should be 2)
                keys.append(k)

            if(len(keys) != 2):         # Can only have two children
                print("ERROR: INCORRECT NUMBER OF CHILDREN!")
                return

            #print("Adding: ", keys)
            #print(self.tree)
            # Children of the input node
            child1 = tn(str(keys[0]), self.data[node][keys[0]])     # Creating child 1
            child2 = tn(str(keys[1]), self.data[node][keys[1]])     # Creating child 2
            self.tree.find(str(node)).extend([child1,child2])       # Find parent and add children

            self.help_tree(keys[0])                                 # Build off child 1
            self.help_tree(keys[1])                                 # Build off child 2

        return

    # ======================================================================
    # Initializes construction of the tree
    def tree_data(self):
        node = self.data['root']
        self.tree = tn(str(node))                   # Set root node

        keys = []                                   # List of keys for connected nodes
        for k in self.data[node]:                   # Get all the keys (should be 3)
            keys.append(k)

        if(len(keys) != 3):                         # Root can only have 3 children
            print("ERROR: INCORRECT NUMBER OF CHILDREN!")
            return

        child1 = tn(str(keys[0]), self.data[node][keys[0]], parent=self.tree)
        child2 = tn(str(keys[1]), self.data[node][keys[1]], parent=self.tree)
        child3 = tn(str(keys[2]), self.data[node][keys[2]], parent=self.tree)
        self.tree.extend([child1, child2, child3])  # Add Children to tree (root)

        self.help_tree(keys[0])                     # Build off child 1
        self.help_tree(keys[1])                     # Build off child 2
        self.help_tree(keys[2])                     # Build off child 3
        #print(self.tree)

    # ======================================================================
    # Writes tree data to output file
    def write_tree(self, filename="./output/tree.tre"):
        self.tree.write(filename)

    # ======================================================================
    # Writes edge data to output file
    def write_edges(self, filename="./output/edges.txt"):
        with open(filename, 'w', newline='') as file:

            for x in self.tree.preorder():
                if(x.length == None): continue                  # Skip root node

                name = x.name
                if(name in self.header):
                    name = self.header[x.name]            # Convert to index

                input =  x.parent.name+'\t'+name+'\t'+str(x.length)+'\n'
                file.write(input)


    # ======================================================================
    # Runs R script to creates pdf of a tree
    def create_tree(self, filename="./output/tree.pdf"):
        input = "Rscript hw3-plot-edges.r ./output/edges.txt hw3-tip-labels.txt " + filename
        os.system(input)

    # ======================================================================
    # Runs R script to creates pdf of a newick tree
    def creat_newick(self, filename="./output/tree-newick.pdf"):
        input = "Rscript hw3-plot-newick.r ./output/tree.tre hw3-tip-labels.txt " + filename
        os.system(input)

    # ======================================================================
    # Runs R script to creates pdf of a newick tree
    def debug(self):
        print(self.data)
        print(self.tree)
