###########################################################################
# This script will segment PDB according to user provided start and end IDs
# 
# Written by Xingcheng Lin, 12/12/2016;
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

def segmentPDB( inputFile, outputFile, startresID, endresID ):

    # Get current working directory
    pwd = os.getcwd();

    infile = open(inputFile, 'r');
    outfile = open(outputFile, 'w');


    # Read in lines from the file;

    lines = [line.strip() for line in infile];

    infile.close();

    length = len(lines);
    
    # Flag for output records;
    flag = 0;

    for i in my_lt_range(0, length, 1):

        line = lines[i].split();
        try:
            line[4];
        except IndexError:
            continue;
        else:
            if (line[0] == "ATOM"):
                # We get in the ATOM sections;
                if (int(line[4]) == startresID):
                    flag = 1;
                # The next residue of the end Residue ID is the end point;
                elif (int(line[4]) == (endresID+1)):
                    flag = 0;
                
                # If record matched, output;
                if (flag == 1):
                    outfile.write(lines[i] + "\n");


    outfile.write("END" + "\n");
    outfile.close();
    return 


############################################################################

if __name__ == "__main__":
    inputFile = sys.argv[1];
    outputFile = sys.argv[2];
    startresID = int(sys.argv[3]);
    endresID = int(sys.argv[4]);

    segmentPDB( inputFile, outputFile, startresID, endresID );
    print "When the voice of the Silent touches my words,"
    print "I know him and therefore know myself."





