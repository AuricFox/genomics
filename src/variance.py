import sys
import numpy as np
import matplotlib.pyplot as plt

class Variance:
    def __init__(self, sequences=[], header=[]):
        '''
        Initializes the variance analysis of a series of sequences
        
        Parameter(s):
            sequences (List[str]): a series of sequences being analyzed
            header (List[str]): a series of heder info for the corresponding sequence

        Output(s):
            None
        '''

        self.sequences = sequences                  # List of taxa sequences
        self.header = header                        # Correspoding markers for sequences

        self.data = {}
        self.variance = []
        self.intersect = []

        self.calc_var()
    
    # ==============================================================================================================
    def calc_var(self):
        '''
        Calculates fraction of the most common base

        Parameter(s):
            None

        Output(s):
            None
        '''

        for i in range(len(self.sequences[0])):                  # Iterates thru the bases in the sequence (columns)
            count = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}         # Initializing base counts

            for j in range(len(self.sequences)):                 # Iterates thru the sequences (rows)
                # Increment base count if it exists, else initialize with one
                count[self.sequences[j][i]] = count.get(self.sequences[j][i], 0) + 1
            
            max_value = max(count['A'], count['C'], count['G'], count['T'])         # Find the max non-gap character (most common base)
            total = sum(k for k in count.values())                                  # Sum of characters
            MCB = [key for key,value in count.items() if value == max_value][0]

            count['MCB'] = MCB
            count['variance'] = max_value/total
            self.data[i] = count
            self.variance.append(count['variance'])

    # ==============================================================================================================
    def get_var(self):
        '''
        Retrieves the variance from each column and returns it as a list used for plotting

        Parameter(s):
            None
        
        Output(s):
            data (list[float]): the variance of each column in a sequence
        '''

        data = []
        for key in self.data:              # Iterate thru each dictionary
            data.append(self.data[key]['variance'])    # Retrieve each variance and add it to the list

        return data

    # ==============================================================================================================
    def plot_data(self, display:bool=False, return_fig:bool=True, save_fig:bool=False, smooth:bool=True):
        '''
        Plots raw or smoothed data

        Parameter(s):
            display (bool, optional): toggels the display of the plotted figure
            return_fig(bool, optional): toggels whether the figure is returned
            save_fig (bool, optional): saves the plotted figure to a pdf
            smooth (bool, optional): toggels between raw (false) and smooth (true) data for plotting
        
        Output(s):
            figure (image) or None: if return_fig is true, a figure of the data is returned
            file (pdf) or None: if save_fig is true the figure will be saved to a pdf
        '''
        
        figure = plt.figure()

        # Plot smooth data
        if smooth:
            y = self.moving_avg()
            x = [i+1 for i in range(len(y))]                # X axis, nucleotide base positions
            plt.title('Raw Data w/Smoothing')
        # Plot raw data
        else:
            y = self.variance
            x = [i+1 for i in range(len(y))]                # X axis, nucleotide base positions
            plt.title('Raw Data w/No Smoothing')

        plt.plot(x, y)
        plt.xlabel('Position in the 16S rRNA gene')
        plt.ylabel('Pct Conserved')
  
        if display: plt.show()                              # Toggle show option
        if return_fig: return figure                        # Toggle return option (Plot figure)
        if save_fig: figure.savefig('./output/plot.pdf')
    
    # ==============================================================================================================
    def moving_avg(self, size:int=60):
        '''
        Coverts the variability to moving averages in order to smooth the data
        NOTE: The default size was selected because the output graph closely matched the one in the documentation
        
        Parameter(s):
            size (int, optional): dictates the number of variances used for calculating the moving average
        
        Output(s):
            mavg (List[float]): calculated moving averages
        '''

        if(size >= len(self.variance)):
            print("ERROR: Window size exceeds array size!")
            print(f'Windox Size: {size}, Array Size: {len(self.variance)}')
            return
        
        y = np.array(self.variance)             # Get list of variances

        mavg = []                               # List of moving averages
        for i in range(len(y) - size):          # Slide thru the list of variances until end (Removes trailing zeros)
            avg = sum(y[i:i+size])/size         # Get average of values within window
            mavg.append(avg)

        return mavg

    # ==============================================================================================================
    def plot_v_regions(self, display:bool=False, return_fig:bool=True, save_fig:bool=False, spacing:int=30, numV:int=6):
        '''
        Plots smoothed data with variance regions
        
        Parameter(s):
            display (bool, optional): toggels the display of the figure
            return_fig (bool, optional): toggels whether the figure is returned
            save_fig (bool, optional): saves the plotted figure to a pdf
            spacing (int, optional): specifies the minimum number of bases needed for a v region (peak)
            numV (int, optional): specifies the number of desired v regions (peaks)
        
        Output(s):
            figure (image) or None: if return_fig is true, a figure of the smoothed data with marked v regions is returned
            file (pdf) or None: if save_fig is true the figure will be saved to a pdf
        '''

        figure = plt.figure()
        y_1 = np.array(self.moving_avg())                               # Variability y
        x = np.array([i+1 for i in range(len(y_1))])                     # X axis, nucleotide base positions

        line = 0.0                                                      # Intersection line used for identifying v regions
        intersect = []                                                  # Intersection points

        while(line <= 1.0):
            y_2 = np.array([line]*len(y_1))                                             # Intersecting y
            intersect = np.argwhere(np.diff(np.sign(y_1 - y_2))).flatten()              # Get intersection points
            intersect = np.insert(intersect, [0,len(intersect)], [0,len(self.variance)])# Add start(0) and end(1,514) points

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
        if save_fig: figure.savefig('./output/v_region_plot.pdf')

    # ==============================================================================================================
    def save_plot(self, plot:bool, smooth:bool, vRegion:bool, spacing:int=30, numV:int=6, filename:str='plot.pdf'):
        '''
        Function used to save plotted figures to a file

        Parameter(s):
            plot (bool): save the raw/smooth plot to a file if true
            smooth (bool): smooth the data if true, else plot the raw data
            vRegion (bool): save the plot with variance regions to a file if true
            spacing (int, optional): specifies the minimum number of bases needed for a v region (peak)
            numV (int, optional): specifies the number of desired v regions (peaks)
            filename (str): name of the file being saved (not the path)

        Output(s):
            A file with the saved plot(s) if variables plot or vRegion are true, else outputs None
        '''

        if plot:
            figure = self.plot_data(
                display=False, 
                return_fig=True, 
                save_fig=False, 
                smooth=smooth
            )

            figure.savefig(f'./output/raw_{filename}')
        
        if vRegion:
            figure = self.plot_v_regions(
                display=False, 
                return_fig=True, 
                save_fig=False, 
                spacing=spacing, 
                numV=numV
            )

            figure.savefig(f'./output/v_region_{filename}')

    # ==============================================================================================================
    def __str__(self):
        '''
        Prints class attributes header and data

        Parameter(s):
            None

        Output(s):
            None
        '''

        string = f'HEADER:\n{self.header}\nData:'

        for key in self.data:
            string += f'\n{self.data[key]}'

        return string


def main():
    return


if __name__ == '__main__':
    main()