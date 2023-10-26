####################################################################################
# This script will help read in the extracted [ pairs ] section the smog.top file and
# remove the interactions in the IDR regions

#
# Written by Xingcheng Lin, 04/06/2017
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


def removeIDR_for_pairs_section(inputFile, outputFile):

    pairs_interactions = np.loadtxt('pairs.dat', comments=';')

    outfile = open(outputFile, 'w')

    histonetailSegment = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887]

    # Remove any pairs that involve parts of the histone tails;
    for i in my_lt_range(0, np.shape(pairs_interactions)[0], 1):
        pair_i = int(pairs_interactions[i][0])
        pair_j = int(pairs_interactions[i][1])
        if (pair_i not in histonetailSegment) and (pair_j not in histonetailSegment):
            outfile.write(str(int(pairs_interactions[i][0])) + "\t" + str(int(pairs_interactions[i][1])) + "\t" + str(int(pairs_interactions[i][2])) + "\t" + str(pairs_interactions[i][3]) + "\t" + str(pairs_interactions[i][4]) + "\t" + str(pairs_interactions[i][5]) + "\t" + str(pairs_interactions[i][6]) + "\n")
    outfile.close()

    return

############################################################################


if __name__ == "__main__":

    inputFile = sys.argv[1]
    outputFile = sys.argv[2]
    removeIDR_for_pairs_section(inputFile, outputFile)

    print("Love is an endless mystery,")
    print("for it has nothing else to explain it.")
