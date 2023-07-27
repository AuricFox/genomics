import matplotlib.pyplot as plt
import numpy as np

'''
   REFERENCE SEQUENCE
  |   | A | C | G | T |
----------------------|
  |   |   |   |   |   |A S
----------------------|L E
A |   |   |   |   |   |I Q
----------------------|G U
C |   |   |   |   |   |N E
----------------------|M N
G |   |   |   |   |   |E C
----------------------|N E
T |   |   |   |   |   |T
-----------------------
'''

class alignment:
    def __init__(self, ref:str, seq:str, gap_pen:int=-2, match_pen:int=-1, ignore:bool=False):
        self.ref = ref                  # Sequence being referenced
        self.seq = seq                  # Sequence 1
        self.gap_pen = gap_pen          # Gap penalty
        self.match_pen = match_pen      # Mismatch penalty
        self.ignore = ignore            # Ignore start and end gaps

        self.results = {
            'score': 0,                 # Alignment score
            'reference': '',            # Alignment of the reference sequence with gaps added
            'sequence': '',             # Alignment of the input sequence with gaps added
            'visual': ''                # Visualization of alignment (for text file)
        }
        
        self.counts = {
            'matches': 0,
            'mismatches': 0,
            'gaps': 0,
            'start gaps': 0,
            'end gaps': 0
        }

        self.s = None                   # Scoring matrix
        self.n = None                   # Neighbor matrix (tracks direction)
        self.build_matrix()
        self.get_alignment()
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Initializes scoring and neighbor matrices
    def build_matrix(self):
        col = len(self.ref) + 1
        row = len(self.seq) + 1

        self.s = np.zeros((row, col))                   # Initialize scoring matrix
        self.n = np.full((row,col), None)               # Initialize neighbor matrix

        # Fill the first row with gap penalties and neighbor information
        for i in range(row):
            self.s[i, 0] = 0 if self.ignore else self.gap_pen * i
            self.n[i, 0] = (i - 1, 0)
        
        # Fill the firstcolumn with gap penalties and neighbor information
        for i in range(col):
            self.s[0, i] = 0 if self.ignore else self.gap_pen * i
            self.n[0, i] = (0, i - 1)

        for n in range(1, row):                         # Iterate thru each element once for alignment
            for m in range(1, col):
                # print("Row: ", n, " Col: ", m)
                self.align_sequence(m, n)               # Computing each score at row x column
    
    # ----------------------------------------------------------------------------------------------------------------------
    """
    Updates the neighbor and scoring matrices with corresponding placement scores and neighbor
    Inputs:
        * col (int): Column index of the current cell.
        * row (int): Row index of the current cell.
    """
    def align_sequence(self, col, row):
        #print("Col: ", col," Row: ", row)

        left = self.s[row][col-1] + self.gap_pen                        # Score from left neighbor                             
        right = self.s[row-1][col] + self.gap_pen                       # Score from right neighbor
        corner = self.s[row-1][col-1]                                   # Score from corner neighbor

        if(self.ref[col-1] == self.seq[row-1]):                         # Match, add point
            corner = corner + 1
        else:                                                           # No match, add penalty
            corner = corner + self.match_pen

        score_index = np.argmax([left, corner, right])                  # Take the best score
        neighbors = [(row, col-1), (row-1, col-1), (row-1, col)]        # [left, corner, right]
        best_neighbor = neighbors[score_index]                          # Get neighbor with the best score

        self.s[row][col] = max(left,corner, right)                      # Add best score to matrix
        self.n[row][col] = best_neighbor                                # Add best neighbor to matrix

    # ----------------------------------------------------------------------------------------------------------------------
    # Universal getter for alignment strings
    def get_alignment(self):

        # Resetting alignment values
        self.results["score"] = 0
        self.results["reference"] = ""
        self.results["sequence"] = ""
        self.results["visual"] = ""

        if(self.ignore): 
            self.get_local_alignment()               # Get alignment that ignores start/end gaps
        else: 
            self.get_global_alignment()              # Get alignment that tracks start/end gaps

        return (self.results)
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Creates alignment strings for sequence reference, sequence 1, and alignment visualization. Also computes total alignment score.
    def get_global_alignment(self):

        if(self.ignore):                            # Matrix must be rebuilt to count start/end gaps
            self.ignore = False
            self.build_matrix()

        row = len(self.seq)
        col = len(self.ref)
        position = self.n[row][col]                 # Start alignment here (n-row, n-col)
        
        while(position != None):
            # print(f"Alignment: {position} Row: {row} Col: {col}")
            self.alignment_string(row, col, position)

            # Incrementing values to next neighbor
            row, col = position[0], position[1]
            position = self.n[row][col]

    # ----------------------------------------------------------------------------------------------------------------------
    # # Creates alignment strings for sequence reference, sequence 1, and alignment visualization while ingoring start/end gaps
    def get_local_alignment(self):

        if not self.ignore:                         # Matrix must be rebuilt to ignore start/end gaps
            self.ignore = True
            self.build_matrix()

        row = len(self.seq)                         # Last Row
        col = np.argmax(self.s[row])                # Get index with max value
        position = self.n[row][col]                 # Get contributing neighbor

        i = len(self.ref)
        while(i > col):                             # Adding end characters to alignment strings
            #print("Index: ", i)
            self.results['reference'] = f"{self.ref[i-1]}{self.results['reference']}"
            self.results['sequence'] = f"_{self.results['sequence']}"
            self.results['visual'] = f" {self.results['visual']}"
            self.counts['end gaps'] += 1                                        # Incrementing start/end gap count
            i -= 1

        while(row != 0 and col != 0):                                           # Adding alignment characters
            #print("Alignment: ", position, " Row: ", row, " Col: ", col, " Score: ",  self.results['score'])
            self.alignment_string(row, col, position)

            # Incrementing values to next neighbor
            row = position[0]
            col = position[1]
            position = self.n[row][col]

        while(position != None):                                                # Adding start characters to alignment strings
            #print("Row: ", row, " Col: ", col, " position: ", position)
            if(col == 0):                                                       # End of reference string reached
                self.results['reference'] = f"_{self.results['reference']}"
                self.results['sequence'] = f"{self.seq[row-1]}{self.results['sequence']}"
            elif(row == 0):                                                     # End of sequence being aligned
                self.results['reference'] = f"{self.ref[col-1]}{self.results['reference']}"
                self.results['sequence'] = f"_{self.results['sequence']}"

            self.results['visual'] = f" {self.results['visual']}"
            self.counts['start gaps'] += 1                                      # Incrementing start/end gap count
            
            # Incrementing values to next neighbor
            row = position[0]
            col = position[1]
            position = self.n[row][col]

    # ----------------------------------------------------------------------------------------------------------------------
    # Appends to alignment strings (self.results['reference'], self.results['sequence'], and self.results['visual'])
    def alignment_string(self, row, col, position):

        # Diagonal alignment (Diagonal neighbor is the best neighbor)
        if(position[0] == row-1 and position[1] == col-1):
            self.results['reference'] = f"{self.ref[col-1]}{self.results['reference']}"
            self.results['sequence'] = f"{self.seq[row-1]}{self.results['sequence']}"

            # Characters from both sequences match (Increment match score)
            if(self.ref[col-1] == self.seq[row-1]):
                self.results['visual'] = f"|{self.results['visual']}"
                self.results['score'] += 1
                self.counts['matches'] += 1                                 # Incrementing match count
            # Characters from both sequences don't match (Add penalty to score)
            else:
                self.results['visual'] = f"X{self.results['visual']}"
                self.results['score'] += self.match_pen
                self.counts['mismatches'] += 1                              # Incrementing mismatch count
            
        # Right alignment (Right neighbor is the best neighbor)
        elif(position[0] == row-1 and position[1] == col):                  # Verticle alignment (sequence 1 gap)
            self.results['reference'] = f"_{self.results['reference']}"
            self.results['sequence'] = f"{self.seq[row-1]}{self.results['sequence']}"
            self.results['visual'] = f" {self.results['visual']}"
            self.results['score'] += self.gap_pen                           # Add gap penalty to score
            self.counts['gaps'] += 1                                        # Incrementing gap count
        # Left alignment (Left neighbor is the best neighbor)
        elif(row == position[0] and position[1] == col-1):                  # Horizontal alignment (sequence 2 gap)
            self.results['reference'] = f"{self.ref[col-1]}{self.results['reference']}"
            self.results['sequence'] = f"_{self.results['sequence']}"
            self.results['visual'] = f" {self.results['visual']}"
            self.results['score'] += self.gap_pen                           # Add gap penalty to score
            self.counts['gaps'] += 1                                        # Incrementing gap count
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Prints score and neighbor matrices (I didn't want to keep typing the prints in main)
    def print_matrices(self):
        print(self.s)
        print(self.n)
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Prints aligned sequences with visiualization (I didn't want to keep typing the prints in main as well)
    '''
    Score: 23
    __C_GATT_TCGGG...
      |  ||| |xx||
    AGCA_ATTCTTCGG...
    '''
    def print_alignment(self):
        print(f"Score: {self.results['score']}")
        print(self.results['reference'])
        print(self.results['visual'])
        print(self.results['sequence'])

    # ----------------------------------------------------------------------------------------------------------------------
    # Writes alignment data to a text file
    def alignment_file(self, kmer):
        filename = f'./output/align/alignment_{kmer}.txt'

        with open(filename, 'w') as f:
            f.write(str( self.results['score']) + '\n')
            f.write(self.results['reference'] + '\n')
            f.write(self.results['visual'] + '\n')
            f.write(self.results['sequence'] + '\n')
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Creates a figure displaying the similarities between the two sequences
    def plot_compare(self, kmer, show_fig=False, save_fig=True):

        plt.clf()                                                           # Clear figure
        count = 0
        for c in range(len(self.results['visual'])):                                      # Plot only matches
            if(self.results['visual'][c] == '|'):
                plt.plot(c, c, '.', color='#570503')
                count += 1

        percent = (count/len(self.seq)) * 100
        plt.title('Comparing Assembled Contig(s) With Spike Protein ({:.2f}% Match)'.format(percent))
        plt.xlabel('Assembled SARS Spike Protein ')
        plt.ylabel('Our Assembled Contig(s)')

        if(show_fig): plt.show()
        if(save_fig): 
            file = './output/align/alignment_{}.png'.format(kmer)
            print("Creating Comparison Plot: ", file)
            plt.savefig(file, dpi=500)

# ==========================================================================================================================
# Testing
def main():
    ref = "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGA"
    seq = "AACCCGCCACCATGTTCGTGTTCCTGGTGCTGCTGCCTCTGGTGTCCAGCCAGTGTGTGAACC"

    #b = alignment(ref, seq, -2, -1, True)
    #b.get_local_alignment()
    #b.get_alignment()
    #b.print_alignment()
    #b.print_matrices()

if __name__ == "__main__":
    main()