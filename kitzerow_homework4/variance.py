# Samuel Kitzerow, kitze012
# Homework 4, Finding variable genomic regions

import sys
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

class Variance:
    def __init__(self, data = [], header = []):
        self.data = data                            # List of taxa sequences
        self.header = header                        # Correspoding markers for sequences
        self.var = []
        self.intersect = []
        self.calc_var(data)
    
    # ==============================================================================================================
    # Calculates fraction of the most common base
    def calc_var(self, data):

        for i in range(len(data[0])):                           # Iterates thru columns
            count = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}         # Initializing base counts

            for j in range(len(data)):                          # Iterates thru rows
                if(data[j][i] in count.keys()):
                    count[data[j][i]] += 1                      # Increment count at specified base

                else:                                           # Key needs to be added to count
                    count[data[j][i]] = 1
            
            max_value = max(count['A'], count['C'], count['G'], count['T'])         # Find the max non-gap character (most common base)
            total = sum(k for k in count.values())                                  # Sum of characters
            MCB = [k for k,v in count.items() if v == max_value][0]

            count['MCB'] = MCB
            count['var'] = max_value/total
            self.var.append(count)

    # ==============================================================================================================
    # Retrieves the variance from each column and returns it as a list
    def get_var(self):
        data = []
        for i in range(len(self.var)):          # Iterate thru each dictionary
            data.append(self.var[i]['var'])     # Retrieve each variance and add it to the list

        return data

    # ==============================================================================================================
    # Plots raw Data
    def plot_raw_data(self, show=True, ret=False):
        figure = plt.figure()
        y = self.get_var()
        x = [i+1 for i in range(len(y))]                    # X axis, nucleotide base positions

        plt.plot(x, y)
  
        plt.xlabel('Position in the 16S rRNA gene')         # X axis name
        plt.ylabel('Pct Conserved')                         # Y axis name
        plt.title('Raw Data w/No Smoothing')
  
        if show: plt.show()                                 # Toggle show option
        if ret: return figure                               # Toggle return option (Plot figure)

    # ==============================================================================================================
    # Plots smooth Data after raw data has been converted
    def plot_smooth_data(self, show=True, ret=False):
        figure = plt.figure()
        y = self.moving_avg()                               # Y axis, get moving averages
        x = [i+1 for i in range(len(y))]                    # X axis, nucleotide base positions

        plt.plot(x, y)
  
        plt.xlabel('Position in the 16S rRNA gene')         # X axis name
        plt.ylabel('Pct Conserved')                         # Y axis name
        plt.title('Raw Data w/Smoothing')
  
        if show: plt.show()                                 # Toggle show option
        if ret: return figure                               # Toggle return option (Plot figure)
    
    # ==============================================================================================================
    # Covert the variability to moving averages in order to smooth the data
    # The default size was selected because the output graph closely matched the one in the hw documentation
    def moving_avg(self, size=60):

        if(size >= len(self.var)):
            print("ERROR: Window size exceeds array size!")
            print("Windox Size: ", size, " Array Size: ", len(self.var))
            return
        
        y = np.array(self.get_var())            # Get list of variances

        mavg = []                               # List of moving averages
        for i in range(len(y) - size):          # Slide thru the list of variances until end (Removes trailing zeros)
            avg = sum(y[i:i+size])/size         # Get average of values within window
            mavg.append(avg)

        return mavg        

    # ==============================================================================================================
    # Function used to write plot figures to a pdf file
    def plot_to_pdf(self, fun, filename='./output/plot.pdf'):

        f = fun(False, True)            # Input function call for plotting smooth and raw data
        f.savefig(filename)             # Write plot to pdf

    # ==============================================================================================================
    # Plots smoothed Data with Intersections
    # The value 0.704 was choosen for gathering clear v regoins without mistaking a small peaks as v regions
    def plot_v_regions(self, show=True, ret=False, limit=0.704):
        figure = plt.figure()
        y_1 = np.array(self.moving_avg())                               # Variability y
        y_2 = np.array([limit]*len(y_1))                                # Intersecting y 
        x = np.array([i+1 for i in range(len(y_1))])                    # X axis, nucleotide base positions

        plt.plot(x, y_1)
        plt.plot(x, y_2)
        self.intersect = np.argwhere(np.diff(np.sign(y_1 - y_2))).flatten()                     # Get intersection points
        self.intersect = np.insert(self.intersect, [0,len(self.intersect)], [0,len(self.var)])  # Add start(0) and end(1,514) points

        plt.xlabel('Position in the 16S rRNA gene')                     # X axis name
        plt.ylabel('Pct Conserved')                                     # Y axis name
        plt.title('Smoothed Data w/Intersection line')
  
        if show: plt.show()                                             # Toggle show option
        if ret: return figure                                           # Toggle return option (Plot figure)

    # ==============================================================================================================
    # Prints class attributes
    def debug(self):
        print("HEADER: ", self.header)
        #print(self.smatrix)
        for i in range(len(self.var)):
            print(self.var[i])