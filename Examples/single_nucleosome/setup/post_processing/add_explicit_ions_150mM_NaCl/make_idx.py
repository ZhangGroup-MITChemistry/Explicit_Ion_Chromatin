###########################################################################
# This script will make the index for each DNA bp
# 
# Written by Xingcheng Lin, 03/20/2023;
###########################################################################

import time;
import subprocess;
import os;
import math;
import sys;
import numpy as np;

################################################
def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step

def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################

def make_idx( outputFile, startAtomID, endAtomID, num_bp_DNA ):

    # Get current working directory
    pwd = os.getcwd();

    outfile = open(outputFile, 'w');

    # The strand and anti-strand of DNA go in the opposite direction;

    for i in my_lt_range(0, num_bp_DNA, 1):
        # Because the DNA base pair starts from 1, but python is 0-base;
        outfile.write("[ DNA_bp" + str(i+1) + " ]" + "\n")
        # if the index is the first DNA bp, the strand does not have the P
        if i == 0:
            outfile.write(str(startAtomID) + " " + str(startAtomID+1) + " " + str(endAtomID) + " " + str(endAtomID-1) + " " + str(endAtomID-2) + "\n")
        # if the index is the last DNA bp, the anti-strand does not have the P
        elif i == 146:
            outfile.write(str(startAtomID+2 + 3*(i-1)) + " " + str(startAtomID+2 + 3*(i-1) + 1) + " " + str(startAtomID+2 + 3*(i-1) + 2) + " " + str(endAtomID - 3*i) + " " + str(endAtomID - 3*i - 1) + "\n")
        else:
            outfile.write(str(startAtomID+2 + 3*(i-1)) + " " + str(startAtomID+2 + 3*(i-1) + 1) + " " + str(startAtomID+2 + 3*(i-1) + 2) + " " + str(endAtomID - 3*i) + " " + str(endAtomID - 3*i - 1) + " " + str(endAtomID - 3*i - 2) + "\n")
            
            
    outfile.close();
    return 


############################################################################

if __name__ == "__main__":
    outputFile = sys.argv[1];
    startAtomID = int(sys.argv[2]);
    endAtomID = int(sys.argv[3]);
    num_bp_DNA = int(sys.argv[4])

    make_idx( outputFile, startAtomID, endAtomID, num_bp_DNA );
    print ("When the voice of the Silent touches my words,")
    print ("I know him and therefore know myself.")

