# Samuel Kitzerow, kitze012
# Homework 3, Phylogeny Inference

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

                row_sum = 0.0
                for x in range(size):                                       # Sum the values in the row
                    row_sum += m[x][j]

                col_sum = 0.0
                for y in range(size):                                       # Sum the values in the column
                    col_sum += m[i][y]

                qmatrix[j][i] = qmatrix[i][j] = (m[i][j])*(size - 2) - row_sum - col_sum    # Calculating Q values

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
        #print("\nA: ", a, ", ", sum_a, " B: ", b, ", ", sum_b, " Diff: ", diff, " S: ", s)

        return (0.5)*(m[a][b]) + (diff / s)                     # Calculating distance from a to u

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
        a, b = loc[0][0], loc[0][1]                         # Index location of first min value
        #print(m, "\n\n", qm)

        d_au = self.get_distance(m, a, b)                   # Calculating distance from a to u
        d_bu = m[a][b] - d_au                               # Calculating distance from b to u
        #print("Min: ", qm[a][b], " d_ab: ", m[a][b], " d_au: ", d_au, " d_bu: ", d_bu)

        u = []
        for i in range(m.shape[0]):
            if(i == a): continue
            d_u = (0.5)*(m[a][i] + m[i][b] - m[a][b])       # Calaculating distances of each taxa to u
            u.append(d_u)

        m = np.delete(m, [a], axis=0)                       # Deleting row a
        m = np.delete(m, [a], axis=1)                       # Deleting column a
        #print("Trimming: \n", m, "\n")

        m[a, 0:m.shape[0]] = u                              # Adding u distances to row (previously b)
        m[0:m.shape[0], a] = u                              # Adding u distances to column (previously b)
        #print("Adding new distances: ", u, "\n", m)

        self.edges[node] = {h[a]: d_au, h[b]: d_bu}         # Branch from a to u
        h[b] = node                                         # Add new node to list
        h.pop(a)                                            # Remove extra label

        return (m, h)                         # Return new dmatrix

    # ----------------------------------------------------------------------------------------------------------
    # Joins all the neighboring taxa together
    # Builds list of branche edges
    def neighbor_joining(self):
        m = self.dmatrix
        h = [i+1 for i in range(len(self.header))]          # Convert taxa lables to numeric values
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

        self.edges[node] = {h[0]:d_vw, h[1]:d_wd, h[2]:d_we}
        #print("d_vw: ", d_vw, " d_wd: ", d_wd, " d_we: ", d_we)
        #print("Edge V: ", edgeV, " Edge D: ", edgeD, " Edge E: ", edgeE)

    # ----------------------------------------------------------------------------------------------------------
    # Prints out values for debugging
    def debug(self):
        print(self.header, "\n\n", self.dmatrix, "\n\n")
        print(self.edges)

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
def write_edges(data, filename="./output/edges.txt"):

    return
# ==============================================================================================================
# Checking q-values 
def check():
    c = (7-2)*(0.5789473684210527)
    row = 0.5789473684210527+0.47368421052631576+0.5263157894736842+0.5789473684210527+0.631578947368421+0.6842105263157895
    col = 0.5789473684210527+0.10526315789473684+0.5263157894736842+0.5263157894736842+0.3684210526315789+0.42105263157894735
    print("Value: ", c - row - col)

# ==============================================================================================================
def main():
    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from fna file
    # python .\phylogeny.py input_file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                          # Extracting fna file name from arguments
        genes = get_data(fna_file)                      # Processing genome data from fna file
        data = Phylogeny(genes[0], genes[1])            # Initializing neighbor joining

    # ----------------------------------------------------------------------------------------------------------
    # Processes sequence data from example fna file
    # python .\phylogeny.py
    elif(len(sys.argv) == 1):
        fna_file = "./example/example.fna"
        genes = get_data(fna_file)
        data = Phylogeny(genes[0], genes[1])

        #data.debug()
        #write_dmatrix(data)
        #print(data.q_matrix(data.dmatrix))
        #stuff = data.join_neighbor(data.dmatrix, data.header, len(data.header) + 1)
        data.neighbor_joining()
        data.debug()
        
    # ----------------------------------------------------------------------------------------------------------
    else:
        print("ERROR: INVALID ARGUMENMTS!")

# ================================================================================================================================================
if __name__ == "__main__":
    main()
    