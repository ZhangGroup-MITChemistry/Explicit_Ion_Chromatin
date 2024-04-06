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

#include <math.h>
#include <stdlib.h>
#include "pair_lj12_cut_coul_debye.h"
#include "atom.h"
#include "neigh_list.h"
#include "force.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include <unistd.h>

#include <string.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJ12CutCoulDebye::PairLJ12CutCoulDebye(LAMMPS *lmp) : PairLJCutCoulCut(lmp) {}

/* ---------------------------------------------------------------------- */

PairLJ12CutCoulDebye::~PairLJ12CutCoulDebye()
{
  if (allocated) {
    memory->destroy(exclusion_list);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJ12CutCoulDebye::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double r,rinv,screening;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  int *tag = atom->tag;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype] && exclusion_list[tag[i]][tag[j]] != 1) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq[itype][jtype]) {
          r = sqrt(rsq);
          rinv = 1.0/r;
          screening = exp(-kappa*r);
          forcecoul = qqrd2e * qtmp*q[j] * screening * (kappa + rinv);
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj = r6inv * lj1[itype][jtype]*r6inv ;
        } else forcelj = 0.0;

        fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

//        if (itype>=1 && itype<=14 && jtype>=1 && jtype<=14){
//          printf("LJ12 DNA-DNA: The interaction between type %d and type %d is: %f\n", itype, jtype, fpair);
//        }else if (itype>=15 && itype<=34 && jtype>=15 && jtype<=34){
//          printf("LJ12 Prot-Prot: The interaction between type %d and type %d is: %f\n", itype, jtype, fpair);
//        }else{
//          //printf("The interaction between type %d and type %d is: %f\n", itype, jtype, fpair);
//          //        }
//          //
//        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype])
            ecoul = factor_coul * qqrd2e * qtmp*q[j] * rinv * screening;
          else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            //evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            evdwl = r6inv*lj3[itype][jtype]*r6inv -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJ12CutCoulDebye::settings(int narg, char **arg)
{
  if (narg < 3 || narg > 4) error->all(FLERR,"Illegal pair_style command");

  kappa = force->numeric(FLERR,arg[1]);
  cut_lj_global = force->numeric(FLERR,arg[2]);
  if (narg == 3) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(FLERR,arg[3]);

  // define the list
  int i,j;
  int na = atom->natoms;
  memory->create(exclusion_list,na+1, na+1,"pair:exclusion");
  for (i = 1; i <= na; i++)
    for (j = i+1; j <= na; j++) {
        exclusion_list[i][j] = 0;
        exclusion_list[j][i] = 0;
    }

  FILE *fp = force->open_potential(arg[0]);
  char line[1024];
  if (fp == NULL)
    error->all(FLERR,"Cannot open pair exclusion list file");

  // now read and parse pair list file for real
  char *ptr;
  tagint id1, id2, eflag;

  while(fgets(line,1024,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    // get atom ids of pair
    id1 = ATOTAGINT(ptr);
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted pair exclusion list file");
    id2 = ATOTAGINT(ptr);

    exclusion_list[id1][id2] = 1;
    exclusion_list[id2][id1] = 1;

  }
  fclose(fp);

  // reset cutoffs that were previously set from data file

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j] == 1) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ12CutCoulDebye::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ12CutCoulDebye::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&kappa,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLJ12CutCoulDebye::single(int i, int j, int itype, int jtype,
                                  double rsq,
                                  double factor_coul, double factor_lj,
                                  double &fforce)
{
  double r2inv,r6inv,r,rinv,screening,forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq[itype][jtype]) {
    r = sqrt(rsq);
    rinv = 1.0/r;
    screening = exp(-kappa*r);
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] *
      screening * (kappa + rinv);
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    forcelj = r6inv * lj1[itype][jtype]*r6inv ;
  } else forcelj = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv * screening;
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    //philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    philj = r6inv*lj3[itype][jtype]*r6inv -
      offset[itype][jtype];
    eng += factor_lj*philj;
  }

  return eng;
}
