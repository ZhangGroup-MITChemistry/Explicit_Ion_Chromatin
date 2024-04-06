/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gauss6/cut/coul/debye,PairGAUSS6CutCoulDebye)

#else

#ifndef LMP_PAIR_GAUSS6_CUT_COUL_DEBYE_H
#define LMP_PAIR_GAUSS6_CUT_COUL_DEBYE_H

#include "pair_gauss6_cut_coul_cut.h"

namespace LAMMPS_NS {

class PairGAUSS6CutCoulDebye : public PairGAUSS6CutCoulCut {
 public:
  PairGAUSS6CutCoulDebye(class LAMMPS *);
  virtual ~PairGAUSS6CutCoulDebye() {}
  virtual void compute(int, int);
  void settings(int, char **);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double kappa;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
