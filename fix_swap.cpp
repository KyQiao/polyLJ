/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Kaiyao QIAO
------------------------------------------------------------------------- */

#include "fix_swap.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSwap::FixSwap(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    local_swap_atom_list(nullptr),
    random_equal(nullptr), random_unequal(nullptr), c_pe(nullptr) {
  if (narg != 7)
    error->all(FLERR, "Illegal fix atom/swap command");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  ncycles = utils::inumeric(FLERR, arg[4], false, lmp);
  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  Tlow = utils::numeric(FLERR, arg[6], false, lmp);

  // double temperature = utils::numeric(FLERR, arg[6], false, lmp);
  // beta = 1.0 / (force->boltz * temperature);

  if (nevery <= 0)
    error->all(FLERR, "Illegal fix atom/swap command");
  if (ncycles < 0)
    error->all(FLERR, "Illegal fix atom/swap command");
  if (seed <= 0)
    error->all(FLERR, "Illegal fix atom/swap command");

  // random number generator, same for all procs

  random_equal = new RanPark(lmp, seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp, seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  nswap_attempts = 0.0;
  nswap_successes = 0.0;

  atom_swap_nmax = 0;
  local_swap_atom_list = nullptr;

  // set comm size needed by this Fix

  if (atom->q_flag)
    comm_forward = 2;
  else
    comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixSwap::~FixSwap() {
  delete random_equal;
  delete random_unequal;
}

/* ---------------------------------------------------------------------- */

int FixSwap::setmask() {
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */
//initialization before a run(optional)
// check some parameter
/* ---------------------------------------------------------------------- */
void FixSwap::init() {
  char *id_pe = (char *)"thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  // check that no swappable atoms are in atom->firstgroup
  // swapping such an atom might not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit))
        flag = 1;

    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);

    if (flagall)
      error->all(FLERR, "Cannot do atom/swap on atoms in atom_modify first group");
  }
}

double FixSwap::ComputeTemp() {
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double scalar = 0;
  double t = 0.0;
  double dof = domain->dimension * atom->natoms;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }

  double tfactor = force->mvv2e / (dof * force->boltz);

  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  scalar *= tfactor;
  return scalar;
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo swaps
------------------------------------------------------------------------- */

void FixSwap::pre_exchange() {
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep)
    return;

  if (domain->triclinic)
    domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic)
    domain->lamda2x(atom->nlocal + atom->nghost);
  if (modify->n_pre_neighbor)
    modify->pre_neighbor();
  // this build is necessary since last step may not build the list
  neighbor->build(1);

  energy_stored = energy_full();

  int nsuccess = 0;

  double temperature = ComputeTemp();
  beta = 1.0 / (force->boltz * temperature);

  // add a parameter to set a swap onset
  if (temperature <= Tlow) {
    update_swap_atoms_list();
    for (int i = 0; i < ncycles; i++)
      nsuccess += attempt_swap();

    nswap_attempts += ncycles;
    nswap_successes += nsuccess;

    // used to update some parameter?
    energy_full();
  }
  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixSwap::attempt_swap() {
  if (nswap == 0)
    return 0;
  double *q = atom->q;
  double energy_before = energy_stored;

  int i = -1, j = -1;
  pick_swap_atom(&i, &j);
  double qi = 0, qj = 0, qi_t = 0, qj_t = 0;

  // start swap
  // MPI_Bcast should be called from a fixed thread for all thread
  if (i >= 0) {
    qi_t = q[i];
    // MPI_Bcast(&qi, 1, MPI_DOUBLE, comm->me, world);
  }
  if (j >= 0) {
    qj_t = q[j];
    // MPI_Bcast(&qj, 1, MPI_DOUBLE, comm->me, world);
  }
  MPI_Allreduce(&qi_t, &qi, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&qj_t, &qj, 1, MPI_DOUBLE, MPI_SUM, world);
  if (i >= 0) {
    if (atom->q_flag)
      atom->q[i] = qj;
  }
  if (j >= 0) {
    if (atom->q_flag)
      atom->q[j] = qi;
  }
  // start communication with other threads
  // give swap info to other thread
  // no need to build the neighbor list since all
  // possible neighbor are included
  comm->forward_comm_fix(this);

  double energy_after = energy_full();

  if (random_equal->uniform() < exp(beta * (energy_before - energy_after))) {
    // success
    // no need to update list
    // no particle will move and no type will change
    // update_swap_atoms_list();
    energy_stored = energy_after;
    return 1;
  } else {
    // fail
    // swap back
    if (i >= 0) {
      if (atom->q_flag)
        atom->q[i] = qi;
    }
    if (j >= 0) {
      if (atom->q_flag)
        atom->q[j] = qj;
    }
    energy_stored = energy_before;
    comm->forward_comm_fix(this);
  }
  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixSwap::energy_full() {
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_neighbor)
    modify->pre_neighbor();
  if (modify->n_pre_force)
    modify->pre_force(vflag);

  if (force->pair)
    force->pair->compute(eflag, vflag);

  if (atom->molecular) {
    if (force->bond)
      force->bond->compute(eflag, vflag);
    if (force->angle)
      force->angle->compute(eflag, vflag);
    if (force->dihedral)
      force->dihedral->compute(eflag, vflag);
    if (force->improper)
      force->improper->compute(eflag, vflag);
  }

  if (force->kspace)
    force->kspace->compute(eflag, vflag);

  if (modify->n_post_force)
    modify->post_force(vflag);
  if (modify->n_end_of_step)
    modify->end_of_step();

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixSwap::pick_swap_atom(int *i, int *j) {
  int iwhichglobal = static_cast<int>(nswap * random_equal->uniform());
  int jwhichglobal = static_cast<int>(nswap * random_equal->uniform());
  while (iwhichglobal == jwhichglobal)
    jwhichglobal = static_cast<int>(nswap * random_equal->uniform());

  // find out if choosen particle on this proc
  if ((iwhichglobal >= nswap_before) && (iwhichglobal < nswap_before + nswap_local)) {
    int iwhichlocal = iwhichglobal - nswap_before;
    *i = local_swap_atom_list[iwhichlocal];
  }

  if ((jwhichglobal >= nswap_before) && (jwhichglobal < nswap_before + nswap_local)) {
    int jwhichlocal = jwhichglobal - nswap_before;
    *j = local_swap_atom_list[jwhichlocal];
  }
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixSwap::update_swap_atoms_list() {
  // the info is local
  int nlocal = atom->nlocal;
  double **x = atom->x;

  if (atom->nmax > atom_swap_nmax) {
    memory->sfree(local_swap_atom_list);
    // max # of owned+ghost in arrays on this proc
    atom_swap_nmax = atom->nmax;
    local_swap_atom_list = (int *)memory->smalloc(atom_swap_nmax * sizeof(int),
                                                  "SWAP:local_swap_atom_list");
  }
  // # of local swappable particle
  nswap_local = 0;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      local_swap_atom_list[nswap_local] = i;
      nswap_local++;
    }
  }

  MPI_Allreduce(&nswap_local, &nswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&nswap_local, &nswap_before, 1, MPI_INT, MPI_SUM, world);
  nswap_before -= nswap_local;
}

/* ---------------------------------------------------------------------- */

int FixSwap::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;

  if (atom->q_flag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
      buf[m++] = q[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixSwap::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;
  last = first + n;

  if (atom->q_flag) {
    for (i = first; i < last; i++) {
      type[i] = static_cast<int>(buf[m++]);
      q[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++)
      type[i] = static_cast<int>(buf[m++]);
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratio
------------------------------------------------------------------------- */

double FixSwap::compute_vector(int n) {
  if (n == 0)
    return nswap_attempts;
  if (n == 1)
    return nswap_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixSwap::memory_usage() {
  double bytes = (double)atom_swap_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixSwap::write_restart(FILE *fp) {
  int n = 0;
  double list[6];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = nswap_attempts;
  list[n++] = nswap_successes;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), n, fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixSwap::restart(char *buf) {
  int n = 0;
  double *list = (double *)buf;

  seed = static_cast<int>(list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int>(list[n++]);
  random_unequal->reset(seed);

  next_reneighbor = (bigint)ubuf(list[n++]).i;

  nswap_attempts = static_cast<int>(list[n++]);
  nswap_successes = static_cast<int>(list[n++]);

  bigint ntimestep_restart = (bigint)ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR, "Must not reset timestep when restarting fix atom/swap");
}
