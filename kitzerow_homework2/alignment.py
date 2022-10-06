import numpy as np

class alignment:
    def __init__(self, ref, seq, gap_pen = -2, match_pen = -1, ign = False):
        self.ref = ref                          # Sequence being referenced
        self.seq = seq                          # Sequence 1
        self.gap_pen = gap_pen                  # Gap penalty
        self.match_pen = match_pen              # Mismatch penalty
        self.ign = ign                          # Ignore start and end gaps

        self.score = 0                          # Alignment score
        self.ref_align = ""                     # Alignment of sequence reference with gaps added
        self.seq_align = ""                     # Alignment of sequence 1 with gaps added
        self.vis = ""                           # Visualization of alignment (for text file)

        self.s = None                           # Scoring matrix
        self.n = None                           # Neighbor matrix (tracks direction)
        self.build_matrix()
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Setter functions for class attributes
    def set_ref(self, ref):
        self.ref = ref
    def set_seq(self, seq):
        self.seq = seq
    def set_gap_pen(self, gap_pen):
        self.gap_pen = gap_pen
    def set_match_pen(self, match_pen):
        self.match_pen = match_pen

    # ----------------------------------------------------------------------------------------------------------------------
    # Initializes scoring and neighbor matrices
    def build_matrix(self):
        col = len(self.ref) + 1
        row = len(self.seq) + 1

        self.s = np.array([[0]*col for i in range(row)])                # Initialize scoring matrix
        self.n = np.array([[None]*col for i in range(row)])             # Initialize neighbor matrix

        k = 0
        for i in range(col):
            k = (self.gap_pen * i)
            if(self.ign): self.s[0, i] = 0                              # Ignore start and end gap penalties by zeroing
            else: self.s[0, i] = k                                      # Adding row gap penalties
            self.n[0, i] = (0,i-1)                                      # Adding row placeholders for neighbor matrix
        
        k = 0
        for i in range(row):
            k = (self.gap_pen * i)
            if(self.ign): self.s[i, 0] = 0                              # Ignore start and end gap penalties by zeroing
            else: self.s[i, 0] = k                                      # Adding column gap penalties
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
        if(self.ref[col-1] == self.seq[row-1]):                         # Check for matches
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
    # Universal getter for alignment strings
    def get_alignment(self):

        # Resetting alignment values
        self.ref_align = ""
        self.seq_align = ""
        self.vis = ""
        self.score = 0

        if(self.ign): self.get_local_alignment()            # Get alignment that ignores start/end gaps
        else: self.get_global_alignment()                   # Get alignment the tracks start/end gaps
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Creates alignment strings for sequence reference, sequence 1, and alignment visualization. Also computes total alignment score.
    def get_global_alignment(self):

        if(self.ign):                                                           # Matrix must be rebuilt to count start/end gaps
            self.ign = False
            self.build_matrix()

        row = len(self.seq)
        col = len(self.ref)
        pos = self.n[row][col]                                                  # Start alignment here
        
        while(pos != None):
            # print("Alignment: ", pos, " Row: ", row, " Col: ", col)
            self.alignment_string(row, col, pos)

            # Incrementing values to next neighbor
            row = pos[0]
            col = pos[1]
            pos = self.n[row][col]
        
        return (self.ref_align, self.vis, self.seq_align, self.score)

    # ----------------------------------------------------------------------------------------------------------------------
    # # Creates alignment strings for sequence reference, sequence 1, and alignment visualization while ingoring start/end gaps
    def get_local_alignment(self):

        if(self.ign == False):              # Matrix must be rebuilt to ignore start/end gaps
            self.ign = True
            self.build_matrix()

        row = len(self.seq)
        col = np.argmax(self.s[row])
        pos = self.n[row][col]

        #TODO: ADD ends to sequence strings

        while(pos != None):
            print("Alignment: ", pos, " Row: ", row, " Col: ", col)
            self.alignment_string(row, col, pos)

            # Incrementing values to next neighbor
            row = pos[0]
            col = pos[1]
            pos = self.n[row][col]

        #TODO: ADD starts to sequence strings

    # ----------------------------------------------------------------------------------------------------------------------
    # Appends to alignment strings (self.ref_align, self.seq_align, and self.vis)
    def alignment_string(self, row, col, pos):

        if(pos[0] == row-1 and pos[1] == col-1 and self.ref[col-1] == self.seq[row-1]):     # Diagonal alignment (match)
            self.ref_align = self.ref[col-1] + self.ref_align
            self.seq_align = self.seq[row-1] + self.seq_align
            self.vis = "|" + self.vis
            self.score += 1
            
        elif(pos[0] == row-1 and pos[1] == col-1 and self.ref[col-1] != self.seq[row-1]):   # Diagonal alignment (mismatch)
            self.ref_align = self.ref[col-1] + self.ref_align
            self.seq_align = self.seq[row-1] + self.seq_align
            self.vis = "X" + self.vis
            self.score += self.match_pen                                                    # Add mismatch penalty to score

        elif(pos[0] == row-1 and pos[1] == col):                                            # Verticle alignment (sequence 1 gap)
            self.ref_align = "_" + self.ref_align
            self.seq_align = self.seq[row-1] + self.seq_align
            self.vis = " " + self.vis
            self.score += self.gap_pen                                                      # Add gap penalty to score

        elif(row == pos[0] and pos[1] == col-1):                                            # Horizontal alignment (sequence 2 gap)
            self.ref_align = self.ref[col-1] + self.ref_align
            self.seq_align = "_" + self.seq_align
            self.vis = " " + self.vis
            self.score += self.gap_pen                                                      # Add gap penalty to score
    
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
        print(self.ref_align)
        print(self.vis)
        print(self.seq_align)

# ==========================================================================================================================
# Testing
if __name__ == "__main__":
    ref = "GTTACGTCCGTAACC"
    seq = "ACGT"

    a = alignment(ref, seq)
    #a.print_matrices()
    a.get_alignment()
    a.print_alignment()
    #a.print_matrices()

    b = alignment(ref, seq, -2, -1, True)
    b.get_local_alignment()
    #b.get_alignment()
    b.print_alignment()
    b.print_matrices()