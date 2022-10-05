import numpy as np

class alignment:
    def __init__(self, seq1, seq2, gap_pen = -2, match_pen = -1):
        self.seq1 = seq1                        # Sequence 1
        self.seq2 = seq2                        # Sequence 2
        self.gap_pen = gap_pen                  # Gap penalty
        self.match_pen = match_pen              # Mismatch penalty

        self.score = 0                          # Alignment score
        self.seq1_align = ""                    # Alignment of sequence 1 with gaps added
        self.seq2_align = ""                    # Alignment of sequence 2 with gaps added
        self.vis = ""                           # Visualization of alignment (for text file)

        self.s = None                           # Scoring matrix
        self.n = None                           # Neighbor matrix (tracks direction)
        self.build_matrix()
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Setter functions for class attributes
    def set_seq1(self, seq1):
        self.seq1 = seq1
    def set_seq2(self, seq2):
        self.seq2 = seq2
    def set_gap_pen(self, gap_pen):
        self.gap_pen = gap_pen
    def set_match_pen(self, match_pen):
        self.match_pen = match_pen

    # ----------------------------------------------------------------------------------------------------------------------
    # Initializes scoring and neighbor matrices
    def build_matrix(self):
        col = len(self.seq1) + 1
        row = len(self.seq2) + 1

        self.s = np.array([[0]*col for i in range(row)])                # Initialize scoring matrix
        self.n = np.array([[None]*col for i in range(row)])             # Initialize neighbor matrix

        k = 0
        for i in range(col):
            k = (self.gap_pen * i)
            self.s[0, i] = k                                            # Adding row gap penalties
            self.n[0, i] = (0,i-1)                                      # Adding row placeholders for neighbor matrix
        
        k = 0
        for i in range(row):
            k = (self.gap_pen * i)
            self.s[i, 0] = k                                            # Adding column gap penalties
            self.n[i, 0] = (i-1,0)                                      # Adding column placeholders for neighbor matrix

        self.n[0, 0] = None

        for n in range(1, row):                                         # Iterate thru each element once for alignment
            for m in range(1, col):
                # print("Row: ", n, " Col: ", m)
                self.align_sequence(m, n)                               # Computing each score at row x column
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Populates neighbor and scoring matrices with corresponding placement scores and neighbors
    def align_sequence(self, col, row):
        # print("M: ", col," N: ", row)

        k = 0
        if(self.seq1[col-1] == self.seq2[row-1]):                       # Check for matches
            k = 1
        else:                                                           # No match, add penalty
            k = self.match_pen


        left = self.s[row][col-1] + self.gap_pen                        # Score from right neighbor
        corner = self.s[row-1][col-1] + k                               # Score from corner neighbor
        right = self.s[row-1][col] + self.gap_pen                       # Score from right neighbor

        score = max(left, corner, right)                                # Take the best score

        if(left == score):                                              # Best path is left neighbor
            self.n[row][col] = (row, col-1)
        elif(corner == score):                                          # Best path is corner neighbor
            self.n[row][col] = (row-1, col-1)
        else:                                                           # Best path is right neighbor
            self.n[row][col] = (row-1, col)

        self.s[row][col] = score                                        # Add score to matrix

    # ----------------------------------------------------------------------------------------------------------------------
    # Creates alignment strings for sequence 1, sequence 2, and alignment visualization. Also computes total alignment score.
    def get_alignment(self):
        row = len(self.seq2)
        col = len(self.seq1)
        pos = self.n[row][col]                                                                  # Start alignment here
        
        while(pos != None):
            # print("Alignment: ", pos, " Row: ", row, " Col: ", col)

            if(pos[0] == row-1 and pos[1] == col-1 and self.seq1[col-1] == self.seq2[row-1]):   # Diagonal alignment (match)
                self.seq1_align = self.seq1[col-1] + self.seq1_align
                self.seq2_align = self.seq2[row-1] + self.seq2_align
                self.vis = "|" + self.vis
                self.score += 1
            
            elif(pos[0] == row-1 and pos[1] == col-1 and self.seq1[col-1] != self.seq2[row-1]): # Diagonal alignment (mismatch)
                self.seq1_align = self.seq1[col-1] + self.seq1_align
                self.seq2_align = self.seq2[row-1] + self.seq2_align
                self.vis = "X" + self.vis
                self.score += self.match_pen                                                    # Add mismatch penalty to score

            elif(pos[0] == row-1 and pos[1] == col):                                            # Verticle alignment (sequence 1 gap)
                self.seq1_align = "_" + self.seq1_align
                self.seq2_align = self.seq2[row-1] + self.seq2_align
                self.vis = " " + self.vis
                self.score += self.gap_pen                                                      # Add gap penalty to score

            elif(row == pos[0] and pos[1] == col-1):                                            # Horizontal alignment (sequence 2 gap)
                self.seq1_align = self.seq1[col-1] + self.seq1_align
                self.seq2_align = "_" + self.seq2_align
                self.vis = " " + self.vis
                self.score += self.gap_pen                                                      # Add gap penalty to score

            # Incrementing values to next neighbor
            row = pos[0]
            col = pos[1]
            pos = self.n[row][col]
        
        return (self.seq1_align, self.vis, self.seq2_align, self.score)
    
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
        print("Score: ", self.score)
        print(self.seq1_align)
        print(self.vis)
        print(self.seq2_align)

# ==========================================================================================================================
# Testing
if __name__ == "__main__":
    seq1 = "ACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTA"
    seq2 = "ACGTCACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTAACCGTA"

    a = alignment(seq1, seq2)
    #a.print_matrices()
    a.get_alignment()
    a.print_alignment()
    #a.print_matrices()