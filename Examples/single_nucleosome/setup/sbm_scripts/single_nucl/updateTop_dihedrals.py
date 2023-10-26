####################################################################################
# This script will help read in the smog.top file and output a new smog.top file 
# based on the updated strength;

#
# Written by Xingcheng Lin, 04/06/2017
####################################################################################

import math;
import subprocess;
import os;
import math;
import numpy as np;
import sys;

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

def updateTop(inputTopFile, editedTopSection, outputTopFile):


    infile1 = open(inputTopFile, 'r');
    infile2 = open(editedTopSection, 'r');

    outfile = open(outputTopFile, 'w');

    lines1 = [line.rstrip() for line in infile1];
    lines2 = [line.rstrip() for line in infile2];

    length1 = len(lines1);
    
    # Flag for whether to update the corresponding section;
    flag = 0;

    for i in my_lt_range(0, length1, 1):
        if (flag == 0 and lines1[i] == ";ai	aj	ak	al	func	 phi0(deg)       Kd              mult"):
            outfile.write(";ai	aj	ak	al	func	 phi0(deg)       Kd              mult" + "\n");
            flag = 1;
        
        if (flag == 0):
            outfile.write(lines1[i] + "\n");
        elif (flag == 1):
            length2 = len(lines2);
            for j in my_lt_range(0, length2, 1):
                outfile.write(lines2[j] + "\n");
            
            # reset flag;
            flag = 2;

        # When it reaches next section, start copying lines;
        if (flag == 2 and lines1[i] == "[ pairs ]"):
            outfile.write("[ pairs ]" + "\n");
            flag = 0;

    return;

############################################################################

if __name__ == "__main__":

    inputTopFile = sys.argv[1]
    editedTopSection = sys.argv[2]
    outputTopFile = sys.argv[3]

    updateTop(inputTopFile, editedTopSection, outputTopFile);


    print("Love is an endless mystery,")
    print("for it has nothing else to explain it.")

