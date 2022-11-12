# Samuel Kitzerow, kitze012
# Homework 3, Phylogeny Inference
# Creates the edge and tree data structures

from skbio.tree import TreeNode as tn

class Ptree:
    def __init__(self, data):
        self.data = data
        self.tree = None
        slef.edges = None

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

            # Children of the input node
            self.tree.find(node).children += [tn(keys[0], self.data[node][keys[0]])]
            self.tree.find(node).children += [tn(keys[1], self.data[node][keys[1]])]

            help_tree(keys[0])          # Build child 1
            help_tree(keys[1])          # Build child 2

        return

    # ======================================================================
    # Initializes construction of the tree
    def tree_data(self):
        node = self.data['root']
        self.tree = tn(node)                   # Set root node

        keys = []                       # List of keys for connected nodes
        for k in self.data[node]:            # Get all the keys (should be 3)
            keys.append(k)

        if(len(keys) != 3):
            print("ERROR: INCORRECT NUMBER OF CHILDREN!")
            return

        tree.children += [tn(keys[0], self.data[node][keys[0]], parent = tree)]
        tree.children += [tn(keys[1], self.data[node][keys[1]], parent = tree)]
        tree.children += [tn(keys[2], self.data[node][keys[2]], parent = tree)]

        help_tree(keys[0])  # Build child 1
        help_tree(keys[1])  # Build child 2
        help_tree(keys[2])  # Build child 3

        return
