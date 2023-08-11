import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')

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

class Alignment:
    def __init__(self, ref:str, seq:str, gap_pen:int=-2, match_point:int=1, match_pen:int=-1, ignore:bool=False):
        '''
        Initializes alignment of the two sequences and the matrices. Uses local alignment if ignore is true, else
        uses global alignment if ignore is false.

        Parameter(s):
            ref (str): the reference sequence being aligned (sequence 1)
            seq (str): the secondary sequence being aligned with the reference sequence (sequence 2)
            gap_pen (int, default=-2): a penalty imposed when there is a gap in the alignment
            match_point (int, default=1): the point(s) added to the score if there is a match in the alignment
            match_pen (int, default=-1): a penalty imposed if there is a mismatch in the alignment
            ignore (bool, default=False): toggles between ingoring end gaps if true and counting end gaps if false
        
        Output(s): None
        '''

        self.ref = ref
        self.seq = seq
        self.gap_pen = gap_pen
        self.match_point = match_point
        self.match_pen = match_pen
        self.ignore = ignore

        self.results = {
            'score': 0,                 # Alignment score
            'reference': '',            # Alignment of the reference sequence with gaps added
            'sequence': '',             # Alignment of the input sequence with gaps added
            'visual': '',               # Visualization of alignment (for text file)
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
    def reset(self):
        '''
        Resets results dictionary, rebuilds the matrices, and aligns sequences
        
        Parameter(s): None
        
        Output(s): None
        '''
        
        self.results = {
            'score': 0, 'reference': '', 'sequence': '', 'visual': '',
            'matches': 0, 'mismatches': 0, 'gaps': 0, 'start gaps': 0, 'end gaps': 0
        }

        self.build_matrix()
        self.get_alignment()

    # ----------------------------------------------------------------------------------------------------------------------
    def build_matrix(self):
        '''
        Initializes scoring and neighbor matrices

        Parameter(s): None
        
        Output(s): None
        '''

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
    def align_sequence(self, col:int, row:int):
        '''
        Updates the neighbor and scoring matrices with corresponding placement scores and neighbor

        Parameter(s):
            col (int): Column index of the current cell.
            row (int): Row index of the current cell.

        Output(s): None
        '''
        #print(f'Col: {col}, Row: {row}')

        left = self.s[row][col-1] + self.gap_pen                    # Score from left neighbor                             
        right = self.s[row-1][col] + self.gap_pen                   # Score from right neighbor
        corner = self.s[row-1][col-1]                               # Score from corner neighbor

        if(self.ref[col-1] == self.seq[row-1]):                     # Match, add point
            corner = corner + self.match_point
        else:                                                       # No match, add penalty
            corner = corner + self.match_pen

        score_index = np.argmax([left, corner, right])              # Take the best score
        neighbors = [(row, col-1), (row-1, col-1), (row-1, col)]    # [left, corner, right]
        best_neighbor = neighbors[score_index]                      # Get neighbor with the best score

        # Add best score to matrix or zero if numbers negative for local alignment
        if self.ignore:
            self.s[row][col] = max(left,corner, right, 0)
        # Add best score to matrix for global alignment
        else:
            self.s[row][col] = max(left,corner, right)
        self.n[row][col] = best_neighbor                            # Add best neighbor to matrix

    # ----------------------------------------------------------------------------------------------------------------------
    def get_alignment(self):
        '''
        Universal method for performing alignment. Performs local alignment if self.ignore is true, else performs global alignment

        Parameter(s): None
        
        Output(s): None
        '''

        if(self.ignore): 
            self.get_local_alignment()               # Get alignment that ignores start/end gaps
        else: 
            self.get_global_alignment()              # Get alignment that tracks start/end gaps

        return (self.results)
    
    # ----------------------------------------------------------------------------------------------------------------------
    def get_global_alignment(self):
        '''
        Creates alignment strings for sequence reference, sequence 1, and alignment visualization. Also computes total 
        alignment score

        Parameter(s): None

        Output(s): None
        '''

        if(self.ignore):                            # Matrix must be rebuilt to count start/end gaps
            self.ignore = False
            self.reset()

        row = len(self.seq)
        col = len(self.ref)
        neighbor = self.n[row][col]                 # Start alignment here (n-row, n-col)
        
        # Iterate until neighbor is row 0 and column 0
        while(neighbor[0] >= 0 and neighbor[1] >= 0):
            #print(f"Alignment: {neighbor} Row: {row} Col: {col}")
            self.alignment_string(row, col, neighbor)

            # Incrementing values to next neighbor
            row, col = neighbor[0], neighbor[1]
            neighbor = self.n[row][col]

    # ----------------------------------------------------------------------------------------------------------------------
    def get_local_alignment(self):
        '''
        Creates alignment strings for reference sequence, sequence 1, and alignment visualization while ingoring start/end gaps

        Parameter(s): None
        
        Output(s): None
        '''

        if not self.ignore:                         # Matrix must be rebuilt to ignore start/end gaps
            self.ignore = True
            self.reset()
        
        row = len(self.seq)                         # Numbre of rows
        col = len(self.ref)                         # Number of columns
        mri = (row, np.argmax(self.s[row]))         # Index of max value in last row
        mci = (np.argmax(self.s[:,col]), col)       # Index of max value in last column
        max_index = [mri, mci]

        max_row_value = self.s[mri[0]][mri[1]]      # Max value in the row
        max_col_value = self.s[mci[0]][mci[1]]      # Max value in the column
        max_score_index = np.argmax([max_row_value, max_col_value]) # Max value index from both

        index = max_index[max_score_index]          # Index of max value in the score matrix
        neighbor = (row, col)                       # Index of next element in path

        # print(f'neighbor: {neighbor}, Start: {point}, Score: {self.s[index[0]][index[1]]}')
        # print(f'Row: {row}, Column: {col}, Start Row: {index[0]}, Start Column: {index[1]}')
        
        # Iterate thru the end gaps in the last row/col of the matrix
        while (row != index[0] or col != index[1]):
            # print(f"Next Point: {point} Row: {row} Col: {col}")
            
            row, col = neighbor[0], neighbor[1]     # Updating values to next point in path

            if index[0] < row:                      # End gaps in the reference sequence
                neighbor = (row-1, col)
                self.results['end gaps'] += 1
            elif index[1] < col:                    # End gaps in the first sequence
                neighbor = (row, col-1)
                self.results['end gaps'] += 1
            else:                                   # No end gaps (first alignment)
                neighbor = self.n[index[0]][index[1]]

            self.alignment_string(row, col, neighbor, True)


        # Iterate until neighbor is row 0 and column 0
        while(neighbor[0] >= 0 and neighbor[1] >= 0):
            # Incrementing values to next neighbor
            row, col = neighbor[0], neighbor[1]
            neighbor = self.n[row][col]

            # Align sequence and add gap penalties
            if (row != 0 and col != 0):
                self.alignment_string(row, col, neighbor, False)
            # End gap reached, don't add gap penalty
            elif ((row != 0 and col == 0) or (row == 0 and col != 0)):
                self.alignment_string(row, col, neighbor, True)
                self.results['start gaps'] += 1

    # ----------------------------------------------------------------------------------------------------------------------
    def alignment_string(self, row:int, col:int, neighbor, ignore:bool=False):
        '''
        Appends to alignment strings (self.results['reference'], self.results['sequence'], and self.results['visual'])

        Parameter(s):
            col (int): Column index of the current cell.
            row (int): Row index of the current cell.
            neighbor (int, int): the neighboring cell that the current cell is pointing too
            ignore (bool, default=False): ignore end gap penalties

        Output(s):
            A string indicating the alignment type (match, mismatch, or gap)
        '''

        # Diagonal alignment (Diagonal neighbor is the best neighbor)
        if(neighbor[0] == row-1 and neighbor[1] == col-1):
            self.results['reference'] = f"{self.ref[col-1]}{self.results['reference']}"
            self.results['sequence'] = f"{self.seq[row-1]}{self.results['sequence']}"

            # Characters from both sequences match (Increment match score)
            if(self.ref[col-1] == self.seq[row-1]):
                self.results['visual'] = f"|{self.results['visual']}"
                self.results['score'] += self.match_point
                self.results['matches'] += 1
                return "match"
            # Characters from both sequences don't match (Add penalty to score)
            else:
                self.results['visual'] = f"X{self.results['visual']}"
                self.results['score'] += self.match_pen
                self.results['mismatches'] += 1
                return "mismatch"
            
        # Right alignment (Add gap to reference sequence)
        elif(neighbor[0] == row-1 and neighbor[1] == col):
            self.results['reference'] = f"_{self.results['reference']}"
            self.results['sequence'] = f"{self.seq[row-1]}{self.results['sequence']}"
            
        # Left alignment (Add gap to sequence)
        elif(row == neighbor[0] and neighbor[1] == col-1):              # Horizontal alignment (sequence 2 gap)
            self.results['reference'] = f"{self.ref[col-1]}{self.results['reference']}"
            self.results['sequence'] = f"_{self.results['sequence']}"
        
        
        if not ignore: self.results['score'] += self.gap_pen            # Add gap penalty to score
        self.results['visual'] = f" {self.results['visual']}"
        self.results['gaps'] += 1                                       # Incrementing gap count
        
        return "gap"

    # ----------------------------------------------------------------------------------------------------------------------
    # Writes alignment data to a text file
    def alignment_file(self, filename:str='./output/alignment.txt'):
        '''
        Parameter(s):
            filename (str, default=./output/alignment.txt): the name of a file where the alignment data will be saved
        
        Output(s):
            A file with the saved alignment data
        '''

        with open(filename, 'w') as f:
            f.write(str( self.results['score']) + '\n')
            f.write(self.results['reference'] + '\n')
            f.write(self.results['visual'] + '\n')
            f.write(self.results['sequence'] + '\n')
    
    # ----------------------------------------------------------------------------------------------------------------------
    def plot_compare(self, filename:str='./temp/alignment_plot.jpg'):
        '''
        Creates a figure displaying the similarities between the two sequences

        Parameter(s):
            filename (str, default=./temp/alignment_plot.jpg): name of the file where the plot is saved if not None, else saves nothing
        
        Output(s):
            returns a plotted figure as a window display if display is true, and ouputs a saved
                file if filename is not None. If both are false then nothing is returned
        '''

        plt.clf()                                                           # Clear figure
        count = 0
        for c in range(len(self.results['visual'])):                        # Plot only matches
            if(self.results['visual'][c] == '|'):
                plt.plot(c, c, '.', color='#570503')
                count += 1

        percent = (count/len(self.seq)) * 100
        plt.title('Comparing Assembled Contig(s) With Spike Protein ({:.2f}% Match)'.format(percent))
        plt.xlabel('Assembled SARS Spike Protein ')
        plt.ylabel('Our Assembled Contig(s)')

        plt.savefig(filename, dpi=500)

    # ----------------------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Combines the aligment score, aligned sequences, and visiualization into one string for printing
        Score: 23
        __C_GATT_TCGGG...
          |  ||| |xx||
        AGCA_ATTCTTCGG...

        Parameter(s): None
        
        Output(s):
            A string containing the alignment score, reference sequence, sequnece 2, and the alignment visualization
        '''
        return (
            f'Score: {self.results["score"]}\n'
            f'{self.results["reference"]}\n'
            f'{self.results["visual"]}\n'
            f'{self.results["sequence"]}'
        )

# ==========================================================================================================================
# Testing
def main():
    
    seq = "THISLINEISATEST"
    ref = "ISALIGNED"

    b = Alignment(ref=ref, seq=seq, gap_pen=-2, match_point=1, match_pen=-1, ignore=True)
    print(b.results)
    #print(b.s)
    #print(b.n)
    b.print_alignment()
    #b.print_matrices()

if __name__ == "__main__":
    main()