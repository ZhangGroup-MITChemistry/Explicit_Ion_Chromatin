/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eric Simon (Cray)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_list_ca.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondListCA::BondListCA(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondListCA::~BondListCA()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(k2);
  }
}

/* ---------------------------------------------------------------------- */

void BondListCA::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,dr2,de_bond;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  int *tag = atom->tag;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type =tag[i1]+ tag[i2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    dr2 = dr*dr;

    // force & energy

    de_bond = 2.0*k2[type]*dr;
    if (r > 0.0) fbond = -de_bond/r;
    else fbond = 0.0;

    if (eflag) ebond = k2[type]*dr2;

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondListCA::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  //memory->create(r0,n+1,"bond:r0");
  //memory->create(k2,n+1,"bond:k2");
  //memory->create(k3,n+1,"bond:k3");
  //memory->create(k4,n+1,"bond:k4");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void BondListCA::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal bond_style command");

  FILE *fp = fopen(arg[0],"r");
  char line[1024];
  if (fp == NULL)
    error->all(FLERR,"Cannot open bond list file");

  int num = 1;
  while(fgets(line,1024,fp)) ++num;
  rewind(fp);
  //int array_size = 4*num+4; // This array size is currently ad hoc
  // This array size is fine is there is only DNA. But if we position DNA
  // after the protein, the size is too small. 
  int array_size = 100*(num+2); // Now set the same as angle list
  memory->create(r0,array_size,"bond:r0");
  memory->create(k2,array_size,"bond:k2");

  char *ptr;
  int idx, id1, id2;

  // Loop through the rest of the lines
  idx = 0;
  while(fgets(line,1024,fp)) { 
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;
    
    // skip comment lines starting with #
    if (*ptr == '#') continue;

    id1 = atoi(ptr);
    ptr = strtok(NULL," \t\n\r\f");

    // The second site
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted bond list file");
    id2 = atoi(ptr);
    
    // Setting the idx in the base array
    // defining idx this way will allow you to find the bond type using the two atom index directly.
    // There is no redundancy since they are always continuous numbers. 
    idx = id1 + id2;
    if ((idx-1) > array_size)
        error->all(FLERR,"Parameter array in bond_list_ca.cpp is too short!");

    // r0
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted bond list/ca file");
    r0[idx] = force->numeric(FLERR,ptr);
    
    //k2
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted bond list/ca file");
    k2[idx] = force->numeric(FLERR,ptr);

  }

  fclose(fp);

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   there are no coeffs to be set, but we need to update setflag and pretend
------------------------------------------------------------------------- */

void BondListCA::coeff(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR, arg[0],atom->nbondtypes,ilo,ihi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
      setflag[i] = 1;
      count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondListCA::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondListCA::write_restart(FILE *fp)
{
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondListCA::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&k2[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondListCA::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,r0[i],k2[i]);
}

/* ---------------------------------------------------------------------- */

double BondListCA::single(int type, double rsq, int i, int j, double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double dr2 = dr*dr;
  double de_bond = 2*k2[type]*dr;
  if (r > 0.0) fforce = -de_bond/r;
  else fforce = 0.0;
  return (k2[type]*dr2);
}
