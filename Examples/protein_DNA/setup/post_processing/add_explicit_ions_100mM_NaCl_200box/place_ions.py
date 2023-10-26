###########################################################################
# This script will place ions on a lattice
# Input: number of ions, output file, the coordinate of the first ion, grid size,
# chainID_sta, resid_sta, num_existing_atoms, atom_type, charge
# Written by Xingcheng Lin, 11/02/2021;
###########################################################################

import time
import subprocess
import os
import math
import sys
import numpy as np

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


def place_ions(num_of_ions, outputFile, sta_coord_x, sta_coord_y, sta_coord_z, grid_size, chainID_sta, resid_sta, num_existing_atoms, atom_type, charge):

    # Get current working directory
    pwd = os.getcwd()

    outfile = open(outputFile, 'w')

    count = 0
    for idx_x in my_lt_range(1, 13, 1):
        for idx_y  in my_lt_range(1, 13, 1):
            for idx_z in my_lt_range(1, 13, 1):
                coord_x = sta_coord_x + idx_x * grid_size
                coord_y = sta_coord_y + idx_y * grid_size
                coord_z = sta_coord_z + idx_z * grid_size

                count += 1
                if count <= num_of_ions:
                    chainID_sta += 1
                    resid_sta += 1
                    num_existing_atoms += 1
                    outfile.write(str(num_existing_atoms) + " " + str(chainID_sta) + " " + str(resid_sta) + " " + str(atom_type) + " " + str(charge) + " " + str(coord_x) + " " + str(coord_y) + " " + str(coord_z) + "\n")

    
    outfile.close()
    return


############################################################################

if __name__ == "__main__":
    num_of_ions = int(sys.argv[1])
    outputFile = sys.argv[2]
    sta_coord_x = float(sys.argv[3])
    sta_coord_y = float(sys.argv[4])
    sta_coord_z = float(sys.argv[5])
    grid_size = float(sys.argv[6])
    chainID_sta = int(sys.argv[7])
    resid_sta = int(sys.argv[8])
    num_existing_atoms = int(sys.argv[9])
    atom_type = int(sys.argv[10])
    charge = float(sys.argv[11])

    place_ions(num_of_ions, outputFile, sta_coord_x, sta_coord_y, sta_coord_z, grid_size, chainID_sta, resid_sta, num_existing_atoms, atom_type, charge)
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
