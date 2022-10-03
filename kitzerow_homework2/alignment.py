import numpy as np

class alignmant:
    def __init__(self, seq1, seq2, gap_pen = -2, match_pen = -1):
        self.seq1 = seq1                        # Sequence 1
        self.seq2 = seq2                        # Sequence 2
        self.gap_pen = gap_pen                  # Gap penalty
        self.match_pen = match_pen              # Mismatch penalty

        self.score = None                       # Alignment score
        self.seq1_align = None                  # Alignment of sequence 1 with gaps added
        self.seq2_align = None                  # Alignment of sequence 2 with gaps added
        self.vis = None                         # Visualization of alignment (for text file)

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
    # Initializes scoring matrix
    def build_matrix(self):
        col = len(self.seq1)
        row = len(self.seq2)

        self.s = np.array([[0]*(col + 1) for i in range(row + 1)])      # Initialize scoring matrix
        self.n = np.array([[None]*col for i in range(row)])             # Initialize neighbor matrix

        k = 0
        for i in range(col):                                            # Adding row gap penalties
            k = (self.gap_pen * i)
            self.s[0, i] = k
        
        k = 0
        for i in range(row):                                            # Adding column gap penalties
            k = (self.gap_pen * i)
            self.s[i, 0] = k
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Populates scoring matrix
    def align_sequence(self, col, row):

        if(len(seq1) <= (col-1) and len(seq2) > (row)):                 # End of row has been reached, move to next row
            self.align_sequence(1, row+1)
        if(len(seq1) <= (col-1) or len(seq2) <= (row-1)):               # End of matrix has been reached, exit loop
            return

        k = 0
        print("M: ", col," N: ", row)
        if(self.seq1[col-1] == self.seq2[row-1]):                       # Check for matches
            k = 1
        else:                                                           # No match, add penalty
            k = self.match_pen


        left = self.s[row][col-1] - self.gap_pen                        # Score from right neighbor
        corner = self.s[row-1][col-1] + k                               # Score from corner neighbor
        right = self.s[row-1][col] - self.gap_pen                       # Score from right neighbor

        score = max(left, corner, right)                                # Take the best score

        if(left == score):                                              # Best path is left neighbor
            self.n[row-1][col-1] = (row, col-1)
        elif(corner == score):                                          # Best path is corner neighbor
            self.n[row-1][col-1] = (row-1, col-1)
        else:                                                           # Best path is right neighbor
            self.n[row-1][col-1] = (row-1, col)

        self.s[row][col] = score                                        # Add score to matrix

        self.align_sequence(col+1, row)                                 # Move to next element in row

# ==========================================================================================================================
# Testing
if __name__ == "__main__":
    seq1 = "ACGTA"
    seq2 = "ACGT"

    a = alignmant(seq1, seq2)
    #print(a.s)
    #print(a.n)
    a.align_sequence(1,1)
    print(a.s)
    print(a.n)