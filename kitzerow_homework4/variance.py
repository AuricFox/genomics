# Samuel Kitzerow, kitze012
# Homework 4, Finding variable genomic regions

import sys
import numpy as np


class Variance:
    def __init__(self, seq = [], header = []):
        self.smatrix = seq                          # List of taxa sequences
        self.header = header                        # Correspoding markers for sequences
        self. var = []

        self.l_to_m(seq)
    
    # ==============================================================================================================
    # Copies the list of lists to a numpy matrix
    def l_to_m(self, data):

        for i in range(len(data)):                              # Iterates thru rows
            count = {'A':0, 'C':0, 'G':0, '-':0}                # Initializing base counts

            for j in range(len(data[0])):                       # Iterates thru columns
                if(data[i][j] in count.keys()):
                    count[data[i][j]] += 1                      # Increment count at specified base

                else:                                           # Key needs to be added to count
                    count[data[i][j]] = 1
            
            max_value = max(count['A'], count['C'], count['G'], count['T'])         # Find the max non-gap character  
            total = sum(k for k in count.values())                                   # Sum of characters
            key = [k for k,v in count.items() if v == max_value][0]

            count['var'] = {key: max_value/total}
            self.var.append(count)
    
    # ==============================================================================================================
    # Prints class attributes
    def debug(self):
        print("HEADER: ", self.header)
        #print(self.smatrix)
        for i in range(len(self.var)):
            print(self.var[i])