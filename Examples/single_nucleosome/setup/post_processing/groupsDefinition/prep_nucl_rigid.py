import numpy as np
import sys
import os

# prepare lammps input rigid body file 

# ----- INPUT PARAMETERS -----

# rigid_n_bp is the total number of nucleosomal DNA bp that we rigidize
# we rigidize the central rigid_n_bp base pairs in the nucleosomal DNA
rigid_n_bp = 73

# set the additional number of bp on both ends of the fiber
n_bp_end_1 = 0
n_bp_end_2 = 0

# set the number of nucleosomes and nrl
n_nucl = 1 # the total number of nucleosomes
nrl = 147 # the nucleosomal repeat length (if n_nucl = 1 then nrl can be any value, should not affect the output)

# ----- FINISH -----

n_ca_per_histone = 974
n_bp_per_nucl = 147

histone_segment = np.arange(1, n_ca_per_histone + 1)
tail_segment = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887], dtype=int)

def get_dna_atom_id_list(n_nucl, nrl, n_bp_end_1, n_bp_end_2):
    # get a list that the atom_id of atoms in the same bp are put in the same sub-list
    output_list = []
    n_res_per_ssdna = (n_nucl - 1)*nrl + n_bp_per_nucl + n_bp_end_1 + n_bp_end_2 # number of residues in each ssDNA
    tot_n_ca = n_nucl*n_ca_per_histone
    atom_id = tot_n_ca + 1
    for i in range(n_res_per_ssdna):
        output_list.append([])
        if i == 0:
            k = 2
        else:
            k = 3
        for j in range(k):
            output_list[i].append(atom_id)
            atom_id += 1
    for i in range(n_res_per_ssdna):
        if i == 0:
            k = 2
        else:
            k = 3
        for j in range(k):
            output_list[n_res_per_ssdna - i - 1].append(atom_id)
            atom_id += 1
    return output_list

rigid_list = []
dna_atom_id_list = get_dna_atom_id_list(n_nucl, nrl, n_bp_end_1, n_bp_end_2)
for i in range(n_nucl):
    rigid_list.append([])
    # add histone globular domain atoms to rigid_list
    histone_segment_i = histone_segment + i*n_ca_per_histone
    tail_segment_i = tail_segment + i*n_ca_per_histone
    for each in histone_segment_i:
        if each not in tail_segment_i:
            rigid_list[i].append(int(each))
    # add nucleosomal DNA atoms to rigid_list
    central_bp_id = n_bp_end_1 + i*nrl + int((n_bp_per_nucl - 1)/2)
    start_bp_id = central_bp_id - int((rigid_n_bp - 1)/2)
    end_bp_id = start_bp_id + rigid_n_bp - 1
    for j in range(start_bp_id, end_bp_id + 1):
        for each in dna_atom_id_list[j]:
            rigid_list[i].append(int(each))

rigid_core_list = []
for i in range(n_nucl):
    rigid_list[i].sort()
    for each in rigid_list[i]:
        rigid_core_list.append(each)

# write group_rigid.txt
with open('group_rigid_innerlayer.txt', 'w') as output_file:
    for i in range(n_nucl):
        rigid_list[i] = [str(int(each)) for each in rigid_list[i]]
        line = 'group nucl%d id %s \n' % (i + 1, ' '.join(rigid_list[i]))
        output_file.write(line)
    rigid_core_list = [str(int(each)) for each in rigid_core_list]
    line = 'group core id %s \n' % ' '.join(rigid_core_list)
    output_file.write(line)


