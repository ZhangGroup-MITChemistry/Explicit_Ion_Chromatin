####################################################################################

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


def add_fakeData(inputfilename, outputfilename):

    # Output the updated data to a new file
    outfile = open(outputfilename, 'w')
    
    data = np.loadtxt(inputfilename)
    dr = 0.02

    for i in my_lt_range(1, 401, 1):
        j = 551 + i
        r = format(dr * i + 12.0, '.3f')

        outfile.write(str(int(j)) + " " + str(r) + " " + str(0.0) + " " + str(0.0) + "\n")

    outfile.close()


    return

############################################################################


if __name__ == "__main__":

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]
    add_fakeData(inputfilename, outputfilename)
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
