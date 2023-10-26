#!/usr/bin/env python

# view the centered 73 bp of core DNA as rigid body, with the non-tail segments of histone

import numpy as np
from numpy import linalg as LA
import sys
sys.path.append('/Users/xl23/bin/lmp_tools/lib')
import lammps_tools as lmp

histoneSegment = np.arange(1, 975)
#print histoneSegment
tailSegment = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887], dtype=int)

num_CA = 974 # number of C-alpha in each histone

topo = lmp.Data()
topo.read_from_file('../data.prot_dna', peptideFlag=1) # What's peptideFlag?
n_atom = len(topo.atoms) # overall number of atoms

atoms = np.zeros((n_atom,4), dtype=int)
for ia, atom in enumerate(topo.atoms):
    atoms[ia,0] = atom['i'] # coarse grained "atom" index
    atoms[ia,1] = atom['mol_i'] # chain index
    atoms[ia,2] = atom['res_i'] # residue index
    atoms[ia,3] = atom['atom_type_i'] # atom type
    # these are just the same order as the columns in data.prot_dna
    # Remember that if atom_type >= 15 and atom_type <= 34, then this CG atom is amino acid

# Calculate the number of nucleosomes (with the data in numpy array atoms)
atoms_histone = atoms[np.logical_and(atoms[:,3] >= 15, atoms[:,3] <= 34)]
atoms_dna = atoms[np.logical_or(atoms[:,3] <= 14, atoms[:,3] >= 35)]
num_nucl = int(len(atoms_histone) / num_CA) # total number of nucleosomes
num_res = max(atoms[:,2]) # total number of residues
num_atom = max(atoms[:,0]) # total number of CG atoms
nrl = int(((num_res - num_CA * num_nucl) / 2 - 147) / (num_nucl - 1))
atoms_dna_1 = atoms_dna[atoms_dna[:,1] == 8 * num_nucl + 1] # first dna chain
atoms_dna_2 = atoms_dna[atoms_dna[:,1] == 8 * num_nucl + 2] # second dna chain

print('len(atoms_histone) =',len(atoms_histone))
print('len(atoms_dna) =',len(atoms_dna))
print('num_nucl =',num_nucl)
print('nrl =',nrl)

#   ----------------------------------------    #
#   selections
#   ----------------------------------------    #

# get alpha carbon coordinates
# now we also consider the whole molecule (including all the chains)
# mask = atoms[:, 1]<=100000
# nucl1Ind = squeeze(atoms[:, 0]) # a collection of all the CG atom index

#   ----------------------------------------    #
#   output
#   ----------------------------------------    #

fh = open('group_rigid_plumed.txt', 'w')

for nucl_index in range(num_nucl):
    # First write histone atoms that are not among tail segment
    fh.write('nucl%d: CENTER ATOMS=' % (nucl_index + 1))
    atoms_histone_this_nucl = atoms_histone[np.logical_and(atoms_histone[:,1] >= 8 * (nucl_index) + 1, atoms_histone[:,1] <= 8 * (nucl_index + 1))]
    for each in atoms_histone_this_nucl:
        if (each[0] % num_CA) not in tailSegment:
            fh.write('%d,' % each[0])
    # then write DNA that should be rigid
    # do not rigidize 37 bp on both ends of 147 bp core DNA
    # Rigidize bp with index 38-110 among 147 bp core DNA
    for each in atoms_dna_1:
        if ((each[2] >= num_CA * num_nucl + nrl * nucl_index + 38) and (each[2] <= num_CA * num_nucl + nrl * nucl_index + 110)):
            fh.write('%d,' % each[0])
    for each in atoms_dna_2:
        if ((each[2] >= num_res - nrl * nucl_index - 109) and (each[2] <= num_res - nrl * nucl_index - 37)):
            fh.write('%d,' % each[0])
    fh.write('\n')
fh.close()

'''
for inucl, atomInd in enumerate([nucl1Ind]): # atomind is the element of [nucl1Ind], while inucl is the index of the element
    fh.write('nucl%d: CENTER  ATOMS='%(inucl+1))
    natom = 0

    for id in range(len(atomInd)):
        if (topo.atoms[atomInd[id]-1]['res_i'] % num_CA) not in tailSegment and :
            fh.write('%d,'%atomInd[id])
            natom+=1
    fh.write('\n')

fh.close()
'''


fh = open('group_rigid.txt', 'w')

for nucl_index in range(num_nucl):
    # First write histone atoms that are not among tail segment
    fh.write('group nucl%d id ' % (nucl_index + 1))
    atoms_histone_this_nucl = atoms_histone[np.logical_and(atoms_histone[:,1] >= 8 * (nucl_index) + 1, atoms_histone[:,1] <= 8 * (nucl_index + 1))]
    for each in atoms_histone_this_nucl:
        if (each[0] % num_CA) not in tailSegment:
            fh.write('%d ' % each[0])
    # then write DNA that should be rigid
    # do not rigidize 37 bp on both ends of 147 bp core DNA
    # Rigidize bp with index 38-110 among 147 bp core DNA
    for each in atoms_dna_1:
        if ((each[2] >= num_CA * num_nucl + nrl * nucl_index + 76) and (each[2] <= num_CA * num_nucl + nrl * nucl_index + 148)):
            fh.write('%d ' % each[0])
    for each in atoms_dna_2:
        if ((each[2] >= num_res - nrl * nucl_index - 147) and (each[2] <= num_res - nrl * nucl_index - 75)):
            fh.write('%d ' % each[0])
    fh.write('\n')


'''
for inucl, atomInd in enumerate([(nucl1Ind)]):
    fh.write('group nucl%d id '%(inucl+1))
    natom = 0

    for id in range(len(atomInd)):
        if (topo.atoms[atomInd[id]-1]['res_i'] % num_CA) not in tailSegment:
            fh.write('%d '%atomInd[id])
            natom+=1
    fh.write('\n')

    print natom
'''

# Output the core region;
fh.write('group core id ')
for nucl_index in range(num_nucl):
    # First write histone atoms that are not among tail segment
    atoms_histone_this_nucl = atoms_histone[np.logical_and(atoms_histone[:,1] >= 8 * (nucl_index) + 1, atoms_histone[:,1] <= 8 * (nucl_index + 1))]
    for each in atoms_histone_this_nucl:
        if (each[0] % num_CA) not in tailSegment:
            fh.write('%d ' % each[0])
    # then write DNA that should be rigid
    # do not rigidize 37 bp on both ends of 147 bp core DNA
    # Rigidize bp with index 38-110 among 147 bp core DNA
    for each in atoms_dna_1:
        if ((each[2] >= num_CA * num_nucl + nrl * nucl_index + 76) and (each[2] <= num_CA * num_nucl + nrl * nucl_index + 148)):
            fh.write('%d ' % each[0])
    for each in atoms_dna_2:
        if ((each[2] >= num_res - nrl * nucl_index - 147) and (each[2] <= num_res - nrl * nucl_index - 75)):
            fh.write('%d ' % each[0])
fh.write('\n')


'''
natom = 0
for inucl, atomInd in enumerate([(nucl1Ind)]):

    for id in range(len(atomInd)):
        # Only output those residues in the ordered segment as defined above;
        if (topo.atoms[atomInd[id]-1]['res_i'] % num_CA) not in tailSegment:
        
            fh.write('%d '%atomInd[id])
            natom+=1

fh.write('\n')
print natom
'''

fh.close()
