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
    def plot_smooth_data(self, show=True, return_fig=False):
        figure = plt.figure()
        y = self.moving_avg()                               # Y axis, get moving averages
        x = [i+1 for i in range(len(y))]                    # X axis, nucleotide base positions

        plt.plot(x, y)
  
        plt.xlabel('Position in the 16S rRNA gene')         # X axis name
        plt.ylabel('Pct Conserved')                         # Y axis name
        plt.title('Raw Data w/Smoothing')
  
        if show: plt.show()                                 # Toggle show option
        if return_fig: return figure                        # Toggle return option (Plot figure)
    
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
    # show: boolean value that toggels the display of the figure
    # return_fig: boolean value that toggels whether the figure is returned
    # spacing: integer that specifies the min number of bases needed for a v region
    # numV: integer that specifies the number of desired v regions
    def plot_v_regions(self, show=True, return_fig=False, spacing=30, numV=6):
        figure = plt.figure()
        y_1 = np.array(self.moving_avg())                               # Variability y
        x = np.array([i+1 for i in range(len(y_1))])                    # X axis, nucleotide base positions

        line = 0.0                                                      # Intersection line used for identifying v regions
        intersect = []                                                  # Intersection points

        while(line <= 1.0):
            y_2 = np.array([line]*len(y_1))                                             # Intersecting y
            intersect = np.argwhere(np.diff(np.sign(y_1 - y_2))).flatten()              # Get intersection points
            intersect = np.insert(intersect, [0,len(intersect)], [0,len(self.var)])     # Add start(0) and end(1,514) points

            if(len(intersect) >= numV*2):                                               # Check for the number of v regions
                check = True

                for i in range(0, numV*2, 2):                                           # Iterate thru intersecting points
                    #print("I: ", intersect[i+1], " J: ", intersect[i], " Total: ", intersect[i+1] - intersect[i])
                    if(intersect[i+1] - intersect[i] < spacing): check = False          # Check spacing between points (fails if too low)
                
                if(check == True): break                                                # Intersecting points looks good

            line += 0.001                                               # Move the intersecting line

        self.intersect = intersect
        plt.plot(x, y_1, color='black')
        #plt.plot(x, y_2)

        for i in range(0, len(intersect), 2):
            x1 = intersect[i]
            x2 = intersect[i+1]
            #print('X1: ', x1, ' X2: ', x2)
            plt.plot([x1,x2],[line,line], color='red')

        plt.xlabel('Position in the 16S rRNA gene')                     # X axis name
        plt.ylabel('Pct Conserved')                                     # Y axis name
        plt.title('Smoothed Data w/Intersection line')
  
        if show: plt.show()                                             # Toggle show option
        if return_fig: return figure                                           # Toggle return option (Plot figure)

    # ==============================================================================================================
    # Prints class attributes
    def debug(self):
        print("HEADER: ", self.header)
        #print(self.smatrix)
        for i in range(len(self.var)):
            print(self.var[i])