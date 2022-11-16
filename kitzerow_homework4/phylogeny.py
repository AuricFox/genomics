# Samuel Kitzerow, kitze012
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
import sys
import tree
import numpy as np
np.set_printoptions(threshold=sys.maxsize, precision=10, linewidth=np.inf)

class Phylogeny:
    def __init__(self, seq = [], header = []):
        self.seq = seq                  # List of taxa sequences
        self.header = header            # Correspoding markers for sequences
        self.dmatrix = None             # Distance matrix
        self.edges = {}                 # Edges of phylogenic tree
        self.tree = None                # Nodes of phylogenic tree

        self.build_Dmatrix()            # Build distance matrix

    # ----------------------------------------------------------------------------------------------------------
    # Calculates % dissimilarity between two sequences
    def gen_distance(self, seq1, seq2):
        match = 0
        mismatch = 0

        if(len(seq1) == 0 or len(seq2) == 0):           # Sequences cannot be null
            print("ERROR: NULL STRING!")
            return
        if(len(seq1) != len(seq2)):
            print("ERROR: lengths are not equal!")      # Lengths of comparing sequences are not equal

        total = len(seq1)                               # Total number of bases
        for (i,j) in zip(seq1,seq2):
            if(i == j):
                match += 1
            else:
                mismatch += 1

        return (mismatch/total)                         # Return % dissimilarity


    # ----------------------------------------------------------------------------------------------------------
    # Builds distance matrix of the input sequences
    # Input must be a nxn matrix
    # Since the matrix is a reflection along the diagonal, only half of it is iterated thru
    # and the values are reflected to their corresponding coordinates [i,j] = [j,i]
    def build_Dmatrix(self):
        size = len(self.seq)                                            # Size of n x n distance matrix
        self.dmatrix = np.array([[0.0]*size for i in range(size)])      # Initialize distance matrix

        for i in range(size):                                           # Sequence 1
            for j in range(i, size):                                    # Sequence 2 (Compare)
                if(i == j): continue
                dis = self.gen_distance(self.seq[i],self.seq[j])        # Getting % dissimilarity between sequences
                #print(i , ", ", j, ", ", dis)
                self.dmatrix[i][j] = self.dmatrix[j][i] = dis           # Setting % dissimilarity in distance matrix

    # ----------------------------------------------------------------------------------------------------------
    # Calculates the Q matrix used for determining which neighbors to join
    # Input must be a nxn matrix
    # Since the matrix is a reflection along the diagonal, only half of it is iterated thru
    # and the values are reflected to their corresponding coordinates [i,j] = [j,i]
    def q_matrix(self, m):
        if(m.shape[0] != m.shape[1]):                                           # Invalid matrix demensions
            print("ERROR: NOT A N X N MATRIX\n", "Matrix: ", m.shape[0], " X ", m.shape[1])
            return

        size = m.shape[1]
        qmatrix = np.array([[0.0]*size for i in range(size)])               # Initialize Q matrix

        for i in range(0, size):                                            # Iterate thru rows
            for j in range(i+1, size):                                      # Iterate thru cols
                if(i == j): continue                                        # Skip diagonal elements

                a_sum = 0.0
                for x in range(size):                                       # Sum the values in a
                    a_sum += m[x][j]

                b_sum = 0.0
                for y in range(size):                                       # Sum the values in b
                    b_sum += m[i][y]

                qmatrix[j][i] = (m[i][j])*(size - 2) - a_sum - b_sum        # Calculating Q values
                qmatrix[i][j] = (m[i][j])*(size - 2) - a_sum - b_sum        # Calculating Q values

        return qmatrix

    # ----------------------------------------------------------------------------------------------------------
    # Calculates distances from each node
    # Takes a dmatrix and two index values
    # Returns calculated distance from a to u
    def get_distance(self, m, a, b):

        if(a >= m.shape[0] or b >= m.shape[0]):
            print("ERROR: INPUT OUT OF RANGE OF MATRIX!")
            return

        sum_a = sum(m[a][q] for q in range(m.shape[0]))         # Summing distances in a
        sum_b = sum(m[b][q] for q in range(m.shape[0]))         # Summing distances in b
        diff = sum_a - sum_b                                    # Computing difference between a and b
        s = 2 * (m.shape[0] - 2)
        distance = (0.5)*(m[a][b]) + (diff / s)                 # Calculating distance from a to u
        #print("A: ", a, ", ", sum_a, " B: ", b, ", ", sum_b, )
        #print("Diff: ", diff, " S: ", s)
        #print("Distance: ", distance, '\n')


        return distance

    # ----------------------------------------------------------------------------------------------------------
    # Joins neighboring sequences
    # Returns a tuple containing (new dmatrix, new list of header info, edge A, edge B)
    def join_neighbor(self, m, h, node):

        if(m.shape[0] != m.shape[1]):
            print("ERROR: MATRIX MUST BE N X N!")
            return

        qm = self.q_matrix(m)                               # Build Q matrix
        min = np.amin(qm)                                   # Min value
        loc = np.where(qm == min)                           # Find mins

        #print('\n', loc)
        #print(m,'\n\n',qm,'\n')

        a, b = loc[0][0], loc[0][1]                         # Index location of first min value
        #print(m, "\n\n", qm)
        #print(h[a], ', ', h[b])
        d_au = self.get_distance(m, a, b)                   # Calculating distance from a to u
        d_bu = m[a][b] - d_au                               # Calculating distance from b to u
        #print("Min: ", qm[a][b], " d_ab: ", m[a][b], " d_au: ", d_au, " d_bu: ", d_bu, '\n')

        u = []
        for i in range(m.shape[0]):
            if(i == a): continue
            d_u = (0.5)*(m[a][i] + m[i][b] - m[a][b])       # Calaculating distances of each taxa to u
            u.append(d_u)

        m = np.delete(m, [a], axis=0)                       # Deleting row a
        m = np.delete(m, [a], axis=1)                       # Deleting column a
        #print("Trimming: \n", m, "\n")

        m[b-1, 0:m.shape[0]] = u                            # Adding u distances to row (b moved 1)
        m[0:m.shape[0], b-1] = u                            # Adding u distances to column (b moved 1)
        #print("Adding new distances: ", u, "\n", m)

        self.edges[node] = {h[a]: d_au, h[b]: d_bu}         # Branch from a to u
        h[b] = node                                         # Add new node to list
        h.pop(a)                                            # Remove extra label

        return (m, h)                         # Return new dmatrix

    # ----------------------------------------------------------------------------------------------------------
    # Joins all the neighboring taxa together
    # Builds list of branche edges
    def neighbor_joining(self):
        m = self.dmatrix.copy()
        h = self.header.copy()                              # taxa lables
        node = len(self.header) + 1                         # Joining node label

        while(m.shape[0] > 3):                              # Iterate thru taxa until there are only two nodes left
            data = self.join_neighbor(m, h, node)           # Join a neighboring taxa
            m = data[0]                                     # Updating matrix
            h = data[1]                                     # Updating header
            node += 1                                       # Updating new node label

            #print(m, '\n', data[2:4], '\n\n')

        #print(self.q_matrix(m))
        d_vw = self.get_distance(m, 0, 1)                   # Distance from v to w
        d_wd = m[0][1] - d_vw                               # Distance from d to w
        d_we = m[0][2] - d_vw                               # Distance from e to w

        self.edges[node] = {h[0]:d_vw, h[1]:d_wd, h[2]:d_we}# Add last three nodes
        self.edges['root'] = node                           # Track root node
        #print("d_vw: ", d_vw, " d_wd: ", d_wd, " d_we: ", d_we)

    # ----------------------------------------------------------------------------------------------------------
    # Prints out values for debugging
    def debug(self):
        #print(self.header, '\n\n', self.dmatrix, '\n\n')
        print(self.edges)

