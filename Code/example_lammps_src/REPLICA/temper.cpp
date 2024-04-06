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
   Contributing author: Mark Sears (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "temper.h"
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "output.h"
#include "thermo.h"
#include "fix.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// #define TEMPER_DEBUG 1

/* ---------------------------------------------------------------------- */

Temper::Temper(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

Temper::~Temper()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_temp;
  delete [] temp2world;
  delete [] world2temp;
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void Temper::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to temper");
  if (domain->box_exist == 0)
    error->all(FLERR,"Temper command before simulation box is defined");

  // lxc editsta: add narg != 8 so that we can deal with two fixes;
  if (narg != 6 && narg != 7 && narg != 8 && narg != 9)
    error->universe_all(FLERR,"Illegal temper command");
  // lxc editend;
  //
  //
  // lxc editsta: add in the option of dealing with 2 fixes;
  // first, if it is in the regular case;
  if (narg == 6 || narg == 7){
  // lxc editend;

    int nsteps = force->inumeric(FLERR,arg[0]);
    nevery = force->inumeric(FLERR,arg[1]);
    double temp = force->numeric(FLERR,arg[2]);
  
    // lxc say: This assigns whichfix to the arg[3], the input we gave for temper fix-ID;
    //
    for (whichfix = 0; whichfix < modify->nfix; whichfix++)
      if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
    // lxc say: If whichfix loops over all the existing fix in the simulation and cannot find
    // one equal to the set fix-ID in the temper command, of course this tempering fix ID is
    // ill-defined
    if (whichfix == modify->nfix)
      error->universe_all(FLERR,"Tempering fix ID is not defined");
  
    seed_swap = force->inumeric(FLERR,arg[4]);
    seed_boltz = force->inumeric(FLERR,arg[5]);
  
    my_set_temp = universe->iworld;
    if (narg == 7) my_set_temp = force->inumeric(FLERR,arg[6]);
  
    // swap frequency must evenly divide total # of timesteps
  
    if (nevery <= 0)
      error->universe_all(FLERR,"Invalid frequency in temper command");
    nswaps = nsteps/nevery;
    if (nswaps*nevery != nsteps)
      error->universe_all(FLERR,"Non integer # of swaps in temper command");
  
    // fix style must be appropriate for temperature control, i.e. it needs
    // to provide a working Fix::reset_target() and must not change the volume.
  
    if ((strcmp(modify->fix[whichfix]->style,"nvt") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/asphere") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/asphere/omp") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/body") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/eff") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/intel") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/kk") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/kk/host") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/kk/device") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/omp") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/sphere") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"nvt/sphere/omp") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"langevin") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"langevin/drude") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"langevin/eff") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"gld") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"gle") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"rigid/nvt") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"rigid/nvt/small") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"rigid/nvt/omp") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"rigid/nvt/small/omp") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/berendsen") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/berendsen/cuda") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/csvr") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/csld") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/rescale") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/rescale/cuda") != 0) &&
        (strcmp(modify->fix[whichfix]->style,"temp/rescale/eff") != 0))
      error->universe_all(FLERR,"Tempering temperature fix is not supported");
  
    // setup for long tempering run
  
    update->whichflag = 1;
    update->nsteps = nsteps;
    update->beginstep = update->firststep = update->ntimestep;
    update->endstep = update->laststep = update->firststep + nsteps;
    if (update->laststep < 0)
      error->all(FLERR,"Too many timesteps");
  
    lmp->init();
  
    // local storage
  
    me_universe = universe->me;
    MPI_Comm_rank(world,&me);
    nworlds = universe->nworlds;
    iworld = universe->iworld;
    boltz = force->boltz;
  
    // pe_compute = ptr to thermo_pe compute
    // notify compute it will be called at first swap
  
    int id = modify->find_compute("thermo_pe");
    if (id < 0) error->all(FLERR,"Tempering could not find thermo_pe compute");
    Compute *pe_compute = modify->compute[id];
    pe_compute->addstep(update->ntimestep + nevery);
  
    // create MPI communicator for root proc from each world
  
    int color;
    if (me == 0) color = 0;
    else color = 1;
    MPI_Comm_split(universe->uworld,color,0,&roots);
  
    // RNGs for swaps and Boltzmann test
    // warm up Boltzmann RNG
  
    if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
    else ranswap = NULL;
    ranboltz = new RanPark(lmp,seed_boltz + me_universe);
    for (int i = 0; i < 100; i++) ranboltz->uniform();
  
    // world2root[i] = global proc that is root proc of world i
  
    world2root = new int[nworlds];
    if (me == 0)
      MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
    MPI_Bcast(world2root,nworlds,MPI_INT,0,world);
  
    // create static list of set temperatures
    // allgather tempering arg "temp" across root procs
    // bcast from each root to other procs in world
  
    set_temp = new double[nworlds];
    if (me == 0) MPI_Allgather(&temp,1,MPI_DOUBLE,set_temp,1,MPI_DOUBLE,roots);
    MPI_Bcast(set_temp,nworlds,MPI_DOUBLE,0,world);
  
    // create world2temp only on root procs from my_set_temp
    // create temp2world on root procs from world2temp,
    //   then bcast to all procs within world
  
    world2temp = new int[nworlds];
    temp2world = new int[nworlds];
    if (me == 0) {
      MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
      for (int i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
    }
    MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);
  
    // if restarting tempering, reset temp target of Fix to current my_set_temp
  
    if (narg == 7) {
      double new_temp = set_temp[my_set_temp];
      modify->fix[whichfix]->reset_target(new_temp);
    }
  
    // setup tempering runs
  
    int i,which,partner,swap,partner_set_temp,partner_world;
    double pe,pe_partner,boltz_factor,new_temp;
  
    if (me_universe == 0 && universe->uscreen)
      fprintf(universe->uscreen,"Setting up tempering ...\n");
  
    update->integrate->setup();
  
    if (me_universe == 0) {
      if (universe->uscreen) {
        fprintf(universe->uscreen,"Step");
        for (int i = 0; i < nworlds; i++)
          fprintf(universe->uscreen," T%d",i);
        fprintf(universe->uscreen,"\n");
      }
      if (universe->ulogfile) {
        fprintf(universe->ulogfile,"Step");
        for (int i = 0; i < nworlds; i++)
          fprintf(universe->ulogfile," T%d",i);
        fprintf(universe->ulogfile,"\n");
      }
      print_status();
    }
  
    timer->init();
    timer->barrier_start();
  
    for (int iswap = 0; iswap < nswaps; iswap++) {
  
      // run for nevery timesteps
  
      update->integrate->run(nevery);
  
      // compute PE
      // notify compute it will be called at next swap
  
      pe = pe_compute->compute_scalar();
      pe_compute->addstep(update->ntimestep + nevery);
  
      // which = which of 2 kinds of swaps to do (0,1)
  
      if (!ranswap) which = iswap % 2;
      else if (ranswap->uniform() < 0.5) which = 0;
      else which = 1;
  
      // partner_set_temp = which set temp I am partnering with for this swap
  
      if (which == 0) {
        if (my_set_temp % 2 == 0) partner_set_temp = my_set_temp + 1;
        else partner_set_temp = my_set_temp - 1;
      } else {
        if (my_set_temp % 2 == 1) partner_set_temp = my_set_temp + 1;
        else partner_set_temp = my_set_temp - 1;
      }
  
      // partner = proc ID to swap with
      // if partner = -1, then I am not a proc that swaps
  
      partner = -1;
      if (me == 0 && partner_set_temp >= 0 && partner_set_temp < nworlds) {
        partner_world = temp2world[partner_set_temp];
        partner = world2root[partner_world];
      }
  
      // swap with a partner, only root procs in each world participate
      // hi proc sends PE to low proc
      // lo proc make Boltzmann decision on whether to swap
      // lo proc communicates decision back to hi proc
  
      swap = 0;
      if (partner != -1) {
        if (me_universe > partner)
          MPI_Send(&pe,1,MPI_DOUBLE,partner,0,universe->uworld);
        else
          MPI_Recv(&pe_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);
  
        if (me_universe < partner) {
          boltz_factor = (pe - pe_partner) *
            (1.0/(boltz*set_temp[my_set_temp]) -
             1.0/(boltz*set_temp[partner_set_temp]));
          if (boltz_factor >= 0.0) swap = 1;
          else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
        }
  
        if (me_universe < partner)
          MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
        else
          MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);
  
  #ifdef TEMPER_DEBUG
        if (me_universe < partner)
          printf("SWAP %d & %d: yes = %d,Ts = %d %d, PEs = %g %g, Bz = %g %g\n",
                 me_universe,partner,swap,my_set_temp,partner_set_temp,
                 pe,pe_partner,boltz_factor,exp(boltz_factor));
  #endif
  
      }
  
      // bcast swap result to other procs in my world
  
      MPI_Bcast(&swap,1,MPI_INT,0,world);
  
      // rescale kinetic energy via velocities if move is accepted
  
      if (swap) scale_velocities(partner_set_temp,my_set_temp);
  
      // if my world swapped, all procs in world reset temp target of Fix
  
      if (swap) {
        new_temp = set_temp[partner_set_temp];
        modify->fix[whichfix]->reset_target(new_temp);
      }
  
      // update my_set_temp and temp2world on every proc
      // root procs update their value if swap took place
      // allgather across root procs
      // bcast within my world
  
      if (swap) my_set_temp = partner_set_temp;
      if (me == 0) {
        MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
        for (i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
      }
      MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);
  
      // print out current swap status
  
      if (me_universe == 0) print_status();
    }
  
    timer->barrier_stop();
  
    update->integrate->cleanup();
  
    Finish finish(lmp);
    finish.end(1);
  
    update->whichflag = 0;
    update->firststep = update->laststep = 0;
    update->beginstep = update->endstep = 0;
  }
  // lxc editsta: else, if there are rigid body fix being defined
  // syntax: temper N M temp fix-ID seed1 seed2 index rigidbodyfix rigidbodyfix-ID
  else if (narg == 8 || narg == 9){

    int nsteps = force->inumeric(FLERR,arg[0]);
    nevery = force->inumeric(FLERR,arg[1]);
    double temp = force->numeric(FLERR,arg[2]);
  
    // lxc say: This assigns whichfix to the arg[3], the input we gave for temper fix-ID;
    // if narg == 8, the rigidbody fix-ID is arg[7], if narg == 9, the rigidbody fix-ID is arg[8];
    //
    for (whichfix1 = 0; whichfix1 < modify->nfix; whichfix1++)
      if (strcmp(arg[3],modify->fix[whichfix1]->id) == 0) break;
    for (whichfix2 = 0; whichfix2 < modify->nfix; whichfix2++)
      if ((narg == 8 && strcmp(arg[7],modify->fix[whichfix2]->id) == 0) || (narg == 9 && strcmp(arg[8],modify->fix[whichfix2]->id) == 0)) break;
    // lxc say: If whichfix loops over all the existing fix in the simulation and cannot find
    // one equal to the set fix-ID in the temper command, of course this tempering fix ID is
    // ill-defined
    if (whichfix1 == modify->nfix || whichfix2 == modify->nfix)
      error->universe_all(FLERR,"Tempering fix ID is not defined");
  
    seed_swap = force->inumeric(FLERR,arg[4]);
    seed_boltz = force->inumeric(FLERR,arg[5]);
  
    my_set_temp = universe->iworld;
    if (narg == 9) my_set_temp = force->inumeric(FLERR,arg[6]);
  
    // swap frequency must evenly divide total # of timesteps
  
    if (nevery <= 0)
      error->universe_all(FLERR,"Invalid frequency in temper command");
    nswaps = nsteps/nevery;
    if (nswaps*nevery != nsteps)
      error->universe_all(FLERR,"Non integer # of swaps in temper command");
  
    // fix style must be appropriate for temperature control, i.e. it needs
    // to provide a working Fix::reset_target() and must not change the volume.
  
    if ((strcmp(modify->fix[whichfix1]->style,"nvt") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/asphere") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/asphere/omp") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/body") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/eff") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/intel") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/kk") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/kk/host") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/kk/device") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/omp") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/sphere") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"nvt/sphere/omp") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"langevin") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"langevin/drude") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"langevin/eff") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"gld") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"gle") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"rigid/nvt") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"rigid/nvt/small") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"rigid/nvt/omp") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"rigid/nvt/small/omp") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/berendsen") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/berendsen/cuda") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/csvr") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/csld") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/rescale") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/rescale/cuda") != 0) &&
        (strcmp(modify->fix[whichfix1]->style,"temp/rescale/eff") != 0))
      error->universe_all(FLERR,"Tempering temperature fix is not supported");
    if ((strcmp(modify->fix[whichfix2]->style,"nvt") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/asphere") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/asphere/omp") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/body") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/eff") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/intel") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/kk") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/kk/host") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/kk/device") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/omp") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/sphere") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"nvt/sphere/omp") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"langevin") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"langevin/drude") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"langevin/eff") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"gld") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"gle") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"rigid/nvt") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"rigid/nvt/small") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"rigid/nvt/omp") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"rigid/nvt/small/omp") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/berendsen") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/berendsen/cuda") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/csvr") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/csld") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/rescale") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/rescale/cuda") != 0) &&
        (strcmp(modify->fix[whichfix2]->style,"temp/rescale/eff") != 0))
      error->universe_all(FLERR,"Tempering temperature fix is not supported");



    // setup for long tempering run
  
    update->whichflag = 1;
    update->nsteps = nsteps;
    update->beginstep = update->firststep = update->ntimestep;
    update->endstep = update->laststep = update->firststep + nsteps;
    if (update->laststep < 0)
      error->all(FLERR,"Too many timesteps");
  
    lmp->init();
  
    // local storage
  
    me_universe = universe->me;
    MPI_Comm_rank(world,&me);
    nworlds = universe->nworlds;
    iworld = universe->iworld;
    boltz = force->boltz;
  
    // pe_compute = ptr to thermo_pe compute
    // notify compute it will be called at first swap
  
    int id = modify->find_compute("thermo_pe");
    if (id < 0) error->all(FLERR,"Tempering could not find thermo_pe compute");
    Compute *pe_compute = modify->compute[id];
    pe_compute->addstep(update->ntimestep + nevery);
  
    // create MPI communicator for root proc from each world
  
    int color;
    if (me == 0) color = 0;
    else color = 1;
    MPI_Comm_split(universe->uworld,color,0,&roots);
  
    // RNGs for swaps and Boltzmann test
    // warm up Boltzmann RNG
  
    if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
    else ranswap = NULL;
    ranboltz = new RanPark(lmp,seed_boltz + me_universe);
    for (int i = 0; i < 100; i++) ranboltz->uniform();
  
    // world2root[i] = global proc that is root proc of world i
  
    world2root = new int[nworlds];
    if (me == 0)
      MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
    MPI_Bcast(world2root,nworlds,MPI_INT,0,world);
  
    // create static list of set temperatures
    // allgather tempering arg "temp" across root procs
    // bcast from each root to other procs in world
  
    set_temp = new double[nworlds];
    if (me == 0) MPI_Allgather(&temp,1,MPI_DOUBLE,set_temp,1,MPI_DOUBLE,roots);
    MPI_Bcast(set_temp,nworlds,MPI_DOUBLE,0,world);
  
    // create world2temp only on root procs from my_set_temp
    // create temp2world on root procs from world2temp,
    //   then bcast to all procs within world
  
    world2temp = new int[nworlds];
    temp2world = new int[nworlds];
    if (me == 0) {
      MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
      for (int i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
    }
    MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);
  
    // if restarting tempering, reset temp target of Fix to current my_set_temp
  
    if (narg == 9) {
      double new_temp = set_temp[my_set_temp];
      modify->fix[whichfix1]->reset_target(new_temp);
      modify->fix[whichfix2]->reset_target(new_temp);
    }
  
    // setup tempering runs
  
    int i,which,partner,swap,partner_set_temp,partner_world;
    double pe,pe_partner,boltz_factor,new_temp;
  
    if (me_universe == 0 && universe->uscreen)
      fprintf(universe->uscreen,"Setting up tempering ...\n");
  
    update->integrate->setup();
  
    if (me_universe == 0) {
      if (universe->uscreen) {
        fprintf(universe->uscreen,"Step");
        for (int i = 0; i < nworlds; i++)
          fprintf(universe->uscreen," T%d",i);
        fprintf(universe->uscreen,"\n");
      }
      if (universe->ulogfile) {
        fprintf(universe->ulogfile,"Step");
        for (int i = 0; i < nworlds; i++)
          fprintf(universe->ulogfile," T%d",i);
        fprintf(universe->ulogfile,"\n");
      }
      print_status();
    }
  
    timer->init();
    timer->barrier_start();
  
    for (int iswap = 0; iswap < nswaps; iswap++) {
  
      // run for nevery timesteps
  
      update->integrate->run(nevery);
  
      // compute PE
      // notify compute it will be called at next swap
  
      pe = pe_compute->compute_scalar();
      pe_compute->addstep(update->ntimestep + nevery);
  
      // which = which of 2 kinds of swaps to do (0,1)
  
      if (!ranswap) which = iswap % 2;
      else if (ranswap->uniform() < 0.5) which = 0;
      else which = 1;
  
      // partner_set_temp = which set temp I am partnering with for this swap
  
      if (which == 0) {
        if (my_set_temp % 2 == 0) partner_set_temp = my_set_temp + 1;
        else partner_set_temp = my_set_temp - 1;
      } else {
        if (my_set_temp % 2 == 1) partner_set_temp = my_set_temp + 1;
        else partner_set_temp = my_set_temp - 1;
      }
  
      // partner = proc ID to swap with
      // if partner = -1, then I am not a proc that swaps
  
      partner = -1;
      if (me == 0 && partner_set_temp >= 0 && partner_set_temp < nworlds) {
        partner_world = temp2world[partner_set_temp];
        partner = world2root[partner_world];
      }
  
      // swap with a partner, only root procs in each world participate
      // hi proc sends PE to low proc
      // lo proc make Boltzmann decision on whether to swap
      // lo proc communicates decision back to hi proc
  
      swap = 0;
      if (partner != -1) {
        if (me_universe > partner)
          MPI_Send(&pe,1,MPI_DOUBLE,partner,0,universe->uworld);
        else
          MPI_Recv(&pe_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);
  
        if (me_universe < partner) {
          boltz_factor = (pe - pe_partner) *
            (1.0/(boltz*set_temp[my_set_temp]) -
             1.0/(boltz*set_temp[partner_set_temp]));
          if (boltz_factor >= 0.0) swap = 1;
          else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
        }
  
        if (me_universe < partner)
          MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
        else
          MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);
  
  #ifdef TEMPER_DEBUG
        if (me_universe < partner)
          printf("SWAP %d & %d: yes = %d,Ts = %d %d, PEs = %g %g, Bz = %g %g\n",
                 me_universe,partner,swap,my_set_temp,partner_set_temp,
                 pe,pe_partner,boltz_factor,exp(boltz_factor));
  #endif
  
      }
  
      // bcast swap result to other procs in my world
  
      MPI_Bcast(&swap,1,MPI_INT,0,world);
  
      // rescale kinetic energy via velocities if move is accepted
  
      if (swap) scale_velocities(partner_set_temp,my_set_temp);
  
      // if my world swapped, all procs in world reset temp target of Fix
  
      if (swap) {
        new_temp = set_temp[partner_set_temp];
        modify->fix[whichfix1]->reset_target(new_temp);
        modify->fix[whichfix2]->reset_target(new_temp);
      }
  
      // update my_set_temp and temp2world on every proc
      // root procs update their value if swap took place
      // allgather across root procs
      // bcast within my world
  
      if (swap) my_set_temp = partner_set_temp;
      if (me == 0) {
        MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
        for (i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
      }
      MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);
  
      // print out current swap status
  
      if (me_universe == 0) print_status();
    }
  
    timer->barrier_stop();
  
    update->integrate->cleanup();
  
    Finish finish(lmp);
    finish.end(1);
  
    update->whichflag = 0;
    update->firststep = update->laststep = 0;
    update->beginstep = update->endstep = 0;
  }
  // lxc editend;
}

/* ----------------------------------------------------------------------
   scale kinetic energy via velocities a la Sugita
------------------------------------------------------------------------- */

void Temper::scale_velocities(int t_partner, int t_me)
{
  double sfactor = sqrt(set_temp[t_partner]/set_temp[t_me]);

  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    v[i][0] = v[i][0]*sfactor;
    v[i][1] = v[i][1]*sfactor;
    v[i][2] = v[i][2]*sfactor;
  }
}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void Temper::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->uscreen," %d",world2temp[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile," %d",world2temp[i]);
    fprintf(universe->ulogfile,"\n");
    fflush(universe->ulogfile);
  }
}
