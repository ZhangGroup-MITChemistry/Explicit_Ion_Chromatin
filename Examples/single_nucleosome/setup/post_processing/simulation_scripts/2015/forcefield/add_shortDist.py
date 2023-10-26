####################################################################################
# This script will help make up for the short distance of the table with
# some artificial energies and forces
# The idea is to take the first row of the table and generate the same value for
# the energies and forces

# Input:  file with a copy of the table
#
# Written by Xingcheng Lin, 09/16/2022
####################################################################################

import math
import subprocess
import os
import math
import numpy as np
import sys

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
###########################################


def add_shortDist(inputfilename, outputfilename):

    data = np.loadtxt(inputfilename)
    data_orig = data[:,1:]
    
    first_energy = data[0, 2]
    first_force = data[0, 3]

    dr = 0.02

    # Add another 50 data points;
    
    new_data = np.empty([49, 3])
    
    for i in my_lt_range(0, 49, 1):
        # We will need to start with 0.02, not 0.0
        idx = i + 1
        # reverse the list
        j = 50 - idx
        r = idx * dr
        diff_r = j * dr
        energy = first_energy + first_force * diff_r

        # The three columns: r, energy, force
        new_data[i, 0] = r
        new_data[i, 1] = energy
        new_data[i, 2] = first_force

    # Concatenate the new data block and the old data block together;

    data_conc = np.vstack((new_data, data_orig))

    # Add back the index to the first column
    index = np.arange(1, np.shape(data_conc)[0]+1).reshape((np.shape(data_conc)[0], 1))

    data_updated = np.concatenate((index, data_conc), axis=1)

    # Output the updated data to a new file
    outfile = open(outputfilename, 'w')

    for i in my_lt_range(0, np.shape(data_updated)[0], 1):
        outfile.write(str(int(data_updated[i, 0])) + " " + str(data_updated[i, 1]) + " " + str(data_updated[i, 2]) + " " + str(data_updated[i, 3]) + "\n")

    outfile.close()


    return

############################################################################


if __name__ == "__main__":

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]
    add_shortDist(inputfilename, outputfilename)
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
