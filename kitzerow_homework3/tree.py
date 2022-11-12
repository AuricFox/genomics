# Samuel Kitzerow, kitze012
# Homework 3, Phylogeny Inference
# Creates the edge and tree data structures

from skbio import TreeNode as tn

class Ptree:
    def __init__(self, data):
        self.data = data
        self.tree = None
        self.edges = None

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

            print(self.tree)
            # Children of the input node
            self.tree.find(str(node)).children += [tn(str(keys[0]), self.data[node][keys[0]])]
            self.tree.find(str(node)).children += [tn(str(keys[1]), self.data[node][keys[1]])]

            self.help_tree(keys[0])          # Build child 1
            self.help_tree(keys[1])          # Build child 2

        return

    # ======================================================================
    # Initializes construction of the tree
    def tree_data(self):
        node = self.data['root']
        self.tree = tn(str(node))                   # Set root node

        keys = []                       # List of keys for connected nodes
        for k in self.data[node]:            # Get all the keys (should be 3)
            keys.append(k)

        if(len(keys) != 3):
            print("ERROR: INCORRECT NUMBER OF CHILDREN!")
            return

        self.tree.children += [tn(str(keys[0]), self.data[node][keys[0]], parent = self.tree)]
        self.tree.children += [tn(str(keys[1]), self.data[node][keys[1]], parent = self.tree)]
        self.tree.children += [tn(str(keys[2]), self.data[node][keys[2]], parent = self.tree)]
        print(self.tree)
        self.help_tree(keys[0])  # Build child 1
        self.help_tree(keys[1])  # Build child 2
        self.help_tree(keys[2])  # Build child 3

        return
