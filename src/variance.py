import sys
import numpy as np
import matplotlib.pyplot as plt

class Variance:
    def __init__(self, sequences = [], header = []):
        self.sequences = sequences                  # List of taxa sequences
        self.header = header                        # Correspoding markers for sequences

        self.data = {}
        self.variance = []
        self.intersect = []

        self.calc_var()
    
    # ==============================================================================================================
    '''
    Calculates fraction of the most common base
    '''
    def calc_var(self):

        for i in enumerate(self.sequences[0]):                  # Iterates thru the bases in the sequence (columns)
            count = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}         # Initializing base counts

            for j in enumerate(self.sequences):                 # Iterates thru the sequences (rows)

                # Increment base count if it exists, else initialize with one
                count[self.sequences[j][i]] = count.get(self.sequences[j][i], 0) + 1
            
            max_value = max(count['A'], count['C'], count['G'], count['T'])         # Find the max non-gap character (most common base)
            total = sum(k for k in count.values())                                  # Sum of characters
            MCB = [key for key,value in count.items() if value == max_value][0]

            count['MCB'] = MCB
            count['variance'] = max_value/total
            self.data[i] = count

    # ==============================================================================================================
    '''
    Retrieves the variance from each column and returns it as a list used for plotting
    Parameter(s): None
    Output(s):
        * data (list[float]): the variance of each column in a sequence
    '''
    def get_var(self):
        data = []
        for key in self.data:              # Iterate thru each dictionary
            data.append(self.data[key]['variance'])    # Retrieve each variance and add it to the list

        return data

    # ==============================================================================================================
    '''
    Plots raw or smoothed data
    Parameter(s):
        * display (bool): toggels the display of the plotted figure
        * return_fig(bool): toggels whether the figure is returned
        * smooth (bool): toggels between raw (false) and smooth (true) data for plotting
    Output(s):
        * figure (image): a figure of the plotted data with marked v regions (peaks)
    '''
    def plot_data(self, display:bool=True, return_fig:bool=False, smooth:bool=True):
        figure = plt.figure()

        # Plot smooth data
        if smooth:
            y = self.moving_avg()                           # Y axis, get moving averages
            x = [i+1 for i in range(len(y))]                # X axis, nucleotide base positions
            plt.title('Raw Data w/Smoothing')
        # Plot raw data
        else:
            y = self.get_var()
            x = [i+1 for i in range(len(y))]                # X axis, nucleotide base positions
            plt.title('Raw Data w/No Smoothing')

        plt.plot(x, y)
        plt.xlabel('Position in the 16S rRNA gene')         # X axis name
        plt.ylabel('Pct Conserved')                         # Y axis name
  
        if display: plt.show()                              # Toggle show option
        if return_fig: return figure                        # Toggle return option (Plot figure)
    
    # ==============================================================================================================
    '''
    Coverts the variability to moving averages in order to smooth the data
    NOTE: The default size was selected because the output graph closely matched the one in the documentation
    Parameter(s):
        * size (int): dictates the number of variances used for calculating the moving average
    Output(s):
        * mavg (List[float]): calculated moving averages
    '''
    def moving_avg(self, size:int=60):

        if(size >= len(self.var)):
            print("ERROR: Window size exceeds array size!")
            print(f'Windox Size: {size}, Array Size: {len(self.var)}')
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
    '''
    Plots smoothed data with variance regions
    Parameter(s):
        * display (bool): toggels the display of the figure
        * return_fig(bool): toggels whether the figure is returned
        * spacing (int): specifies the minimum number of bases needed for a v region (peak)
        * numV (int): specifies the number of desired v regions (peaks)
    Output(s):
        * figure (image): a figure of the smoothed data with marked v regions (peaks)
    '''
    def plot_v_regions(self, display:bool=True, return_fig:bool=False, spacing:int=30, numV:int=6):
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
                    #print(f'I: {intersect[i+1]}, J: {intersect[i]}, Total: {intersect[i+1] - intersect[i]}')
                    if(intersect[i+1] - intersect[i] < spacing): check = False          # Check spacing between points (fails if too low)
                
                if(check == True): break                                                # Intersecting points looks good

            line += 0.001                                               # Move the intersecting line

        self.intersect = intersect
        plt.plot(x, y_1, color='black')
        #plt.plot(x, y_2)

        for i in range(0, len(intersect), 2):
            x1 = intersect[i]
            x2 = intersect[i+1]
            #print(f'X1: {x1}, X2: {x2}')
            plt.plot([x1,x2],[line,line], color='red')

        plt.xlabel('Position in the 16S rRNA gene')                     # X axis name
        plt.ylabel('Pct Conserved')                                     # Y axis name
        plt.title('Smoothed Data w/Intersection line')
  
        if display: plt.show()                                          # Toggle show option
        if return_fig: return figure                                    # Toggle return option (Plot figure)

    # ==============================================================================================================
    # Prints class attributes
    def __str__(self):
        return f'HEADER: {self.header}, Variance(s):\n{self.var}'


def main():
    return


if __name__ == '__main__':
    main()