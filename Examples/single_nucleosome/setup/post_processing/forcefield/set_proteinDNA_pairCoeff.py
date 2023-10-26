#!/usr/bin/env python
#
# protein-DNA 
#   lj/cut/coul/debye args = kappa cutoff (cutoff2)
#   kappa = inverse of the Debye length (inverse distance units)
#   cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
#   cutoff2 = global cutoff for Coulombic (optional) (distance units)

from numpy import *

#   ---         Determine the debye screening length        ---     #
_KB_ = 1.3806505E-23
_PV_ = 8.8541878176E-22
_NA_ = 6.0221415E23
_EC_ = 1.60217653E-19

temp = 300.0
salt = 150.0

t_target = temp
salt_conc = salt / 1000.0

dielectric = 249.4 - 7.88E-01 * t_target + 7.20E-04 * t_target * t_target;
dielectric *= (1.000 - (0.2551 * salt_conc) + (0.05151 * salt_conc * salt_conc) - (0.006889 * salt_conc * salt_conc * salt_conc));

#print dielectric
dielectric = 78.0

ldby = sqrt(dielectric * _PV_ * _KB_ * temp * 1.0E27 / (2.0 * _NA_ * _EC_ * _EC_ * salt_conc))

print ldby

#   ---         Write the pair coefficients         --- #
#
#   --------------------------------------------------- #
filein = open('proteinDna_pairCoeff.in', 'w')

#filein.write('pair_style hybrid vexcluded 2 3.5 3.5 3spn2 ${T} ${salt} 18.0 lj/cut/coul/debye 1.5 2.5 5.0\n')

# Dna pair; calculated using the arithmetic rule
dnaAtomType = 14
dnaEps   = ones(dnaAtomType) *  0.239006;
dnaSigma = zeros(dnaAtomType, dtype=float)
dnaSigma[1  - 1] =  2.0 * 2.25
dnaSigma[2  - 1] =  2.0 * 3.10
dnaSigma[3  - 1] =  2.0 * 2.70
dnaSigma[4  - 1] =  2.0 * 3.55
dnaSigma[5  - 1] =  2.0 * 2.45
dnaSigma[6  - 1] =  2.0 * 3.20
dnaSigma[7  - 1] =  2.0 * 2.70
dnaSigma[8  - 1] =  2.0 * 3.55
dnaSigma[9  - 1] =  2.0 * 2.45
dnaSigma[10 - 1] =  2.0 * 3.20
dnaSigma[11 - 1] =  2.0 * 2.70
dnaSigma[12 - 1] =  2.0 * 3.55
dnaSigma[13 - 1] =  2.0 * 2.45
dnaSigma[14 - 1] =  2.0 * 3.20

# formula from Lammps manual
# http://lammps.sandia.gov/doc/pair_modify.html 
#   mix arithmetic
#   epsilon_ij = sqrt(epsilon_i * epsilon_j)
#   sigma_ij = (sigma_i + sigma_j) / 2 

filein.write('\n')
filein.write('# dna pair-wise interaction\n')
for idat in range(dnaAtomType):
    for jdat in range(idat, dnaAtomType):
    #for jdat in range(dnaAtomType):
        filein.write('pair_coeff %6d %6d 3spn2 %15.8f %15.8f\n'%(idat+1, jdat+1, \
                                                                 sqrt(dnaEps[idat]*dnaEps[jdat]), \
                                                                 (dnaSigma[idat] + dnaSigma[jdat])/2.0 \
                                                                 ) )
# define protein pair
proteinAtomType = 20
proteinEps   = ones(proteinAtomType) *  0.239006;
proteinSigma = ones(proteinAtomType, dtype=float) * 4.5


# protein DNA pair
# pair_coeff for lj/cut specify the following four parameters, and the last two are optional
    #epsilon (energy units)
    #sigma (distance units)
    #cutoff1 (distance units)
    #cutoff2 (distance units) 

filein.write('\n')
filein.write('# protein dna pair-wise interaction\n')
for dat in range(dnaAtomType):
    for pat in range(proteinAtomType):
        eps     = 0.25 *0.5 /4.184        
        sigma   = 5.7   
        ljcut   = sigma*2.5
        filein.write('pair_coeff %6d %6d lj/cut/coul/debye %15.8f %15.8f %15.8f %15.8f\n'%(dat+1, pat+1+dnaAtomType, \
                                                                 eps, sigma, \
                                                                 ljcut, 40) )
# protein - protein electrostatics
# I don't think we need to explicitly say it is zero; that is the default

# ----------------------------------------- #
# protein lj12 electrostatics

filein.write('\n')
filein.write('# protein protein pair-wise interaction\n')

# protein - protein electrostatics
for pati in range(proteinAtomType):
    for patj in range(pati, proteinAtomType):
        eps     = 0.0597515     # we might want to scale this number by a factor of 2 as we did for all other energies
        sigma   = 4.0
        ljcut   = 15.0
        filein.write('pair_coeff %6d %6d lj12/cut/coul/debye %15.8f %15.8f %15.8f %15.8f\n'%(pati+1+dnaAtomType, patj+1+dnaAtomType, \
                                                                 eps, sigma, \
                                                                 ljcut, 40) )

filein.close()

