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
   Contributing author: Kaiyao QIAO
------------------------------------------------------------------------- */

#include "pair_lj_poly.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJPoly::PairLJPoly(LAMMPS *lmp) : Pair(lmp)
{
  // disable respa calculation
  respa_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJPoly::~PairLJPoly()
{
  if (allocated)
  {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJPoly::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;   //get charge for each atom
  double sigma_ij, qtmp; //sigma_ij used in force calculation,qtmp for store charge
  double sigma_ij2, sigma_ij6, sigma_ij12;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++)
  {
    i = ilist[ii];
    qtmp = q[i]; //get charge
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++)
    {
      j = jlist[jj];

      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      sigma_ij = 0.5 * (q[j] + qtmp) * (1 - non_addtive_para * std::abs(q[j] - qtmp));
      sigma_ij2 = sigma_ij * sigma_ij;
      // here cut will be a dimensionless number which represents cut*sigma_ij
      if (rsq < cutsq[itype][jtype] * sigma_ij2)
      {
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;

        // add sigma ij part here
        sigma_ij6 = sigma_ij2 * sigma_ij2 * sigma_ij2;
        sigma_ij12 = sigma_ij6 * sigma_ij6;

        forcelj = r6inv * (lj1[itype][jtype] * sigma_ij12 * r6inv - lj2[itype][jtype] * sigma_ij6);
        fpair = factor_lj * forcelj * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal)
        {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag)
        {
          double ratio = 1 / cut[itype][jtype];
          double ratio6 = pow(ratio, 6.0);
          double ratio12 = ratio6 * ratio6;
          evdwl = r6inv * (lj3[itype][jtype] * sigma_ij12 * r6inv - lj4[itype][jtype] * sigma_ij6) -
                  // offset[itype][jtype];
                  // offset should be calculate everytime
                  4.0 * epsilon[itype][jtype] * (ratio12 - ratio6);
          evdwl *= factor_lj;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair,
                   evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr)
    virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJPoly::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(epsilon, n + 1, n + 1, "pair:epsilon");
  memory->create(lj1, n + 1, n + 1, "pair:lj1");
  memory->create(lj2, n + 1, n + 1, "pair:lj2");
  memory->create(lj3, n + 1, n + 1, "pair:lj3");
  memory->create(lj4, n + 1, n + 1, "pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJPoly::settings(int narg, char **arg)
{
  if (narg != 1)
    error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated)
  {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJPoly::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  // the non_addtive_para can only be set once
  if(comm->me==0){
  if (non_addtive_para_flag == 0)
    non_addtive_para = utils::numeric(FLERR, arg[3], false, lmp);
  non_addtive_para_flag = 1;
  }
  MPI_Bcast(&non_addtive_para_flag, 1, MPI_INT,0, world);
  MPI_Bcast(&non_addtive_para, 1, MPI_DOUBLE, 0, world);
  double cut_one = cut_global;
  if (narg == 5)
    cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++)
  {
    for (int j = MAX(jlo, i); j <= jhi; j++)
    {
      epsilon[i][j] = epsilon_one;
      // sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJPoly::init_style()
{
  // request regular or rRESPA neighbor list

  int irequest;
  int respa = 0;

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style, "^respa"))
  {
    if (((Respa *)update->integrate)->level_inner >= 0)
      respa = 1;
    if (((Respa *)update->integrate)->level_middle >= 0)
      respa = 2;
  }

  irequest = neighbor->request(this, instance_me);

  if (respa >= 1)
  {
    neighbor->requests[irequest]->respaouter = 1;
    neighbor->requests[irequest]->respainner = 1;
  }
  if (respa == 2)
    neighbor->requests[irequest]->respamiddle = 1;

  // set rRESPA cutoffs
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJPoly::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
  {
    // basically this should be set manully
    epsilon[i][j] = (epsilon[i][i] + epsilon[j][j]) / 2.0;
    // the dimensionless cutoff should be the same for the whole system
    cut[i][j] = (cut[i][i] + cut[j][j]) / 2.0;
  }

  // delete sigma part, move that into function compute
  lj1[i][j] = 48.0 * epsilon[i][j];
  lj2[i][j] = 24.0 * epsilon[i][j];
  lj3[i][j] = 4.0 * epsilon[i][j];
  lj4[i][j] = 4.0 * epsilon[i][j];

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce
  // the tail_flag is used only npt or nph is used.
  // for poly disperse system this is not correct
  if (tail_flag)
  {
    // int *type = atom->type;
    // int nlocal = atom->nlocal;

    // double count[2], all[2];
    // count[0] = count[1] = 0.0;
    // for (int k = 0; k < nlocal; k++)
    // {
    //   if (type[k] == i)
    //     count[0] += 1.0;
    //   if (type[k] == j)
    //     count[1] += 1.0;
    // }
    // MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    // double sig2 = sigma[i][j] * sigma[i][j];
    // double sig6 = sig2 * sig2 * sig2;
    // double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
    // double rc6 = rc3 * rc3;
    // double rc9 = rc3 * rc6;
    // etail_ij = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] *
    //            sig6 * (sig6 - 3.0 * rc6) / (9.0 * rc9);
    // ptail_ij = 16.0 * MY_PI * all[0] * all[1] * epsilon[i][j] *
    //            sig6 * (2.0 * sig6 - 3.0 * rc6) / (9.0 * rc9);
    throw std::runtime_error("long-range tail correction is not implemented for polydisperse system");
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJPoly::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
    {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j])
      {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJPoly::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
    {
      if (me == 0)
        utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j])
      {
        if (me == 0)
        {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJPoly::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJPoly::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0)
  {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJPoly::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g\n", i, epsilon[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJPoly::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g\n", i, j, epsilon[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJPoly::single(int i, int j, int itype, int jtype, double rsq,
                          double /*factor_coul*/, double factor_lj,
                          double &fforce)
{
  double r2inv, r6inv, forcelj, philj;
  double sigma_ij = 0.5 * (atom->q[j] + atom->q[i]) * (1 - non_addtive_para * std::abs(atom->q[j] - atom->q[i]));
  double sigma_ij2 = sigma_ij * sigma_ij;
  double sigma_ij6 = sigma_ij2 * sigma_ij2 * sigma_ij2;
  double sigma_ij12 = sigma_ij6 * sigma_ij6;

  r2inv = 1.0 / rsq;
  r6inv = r2inv * r2inv * r2inv;
  forcelj = r6inv * (lj1[itype][jtype] * sigma_ij12 * r6inv - lj2[itype][jtype] * sigma_ij6);
  fforce = factor_lj * forcelj * r2inv;

  // double ratio = sigma_ij / cut[itype][jtype];
  // cutoff set to be a dimensionless one, so that
  double ratio = 1 / cut[itype][jtype];
  double ratio6 = pow(ratio, 6.0);
  double ratio12 = ratio6 * ratio6;
  philj = r6inv * (lj3[itype][jtype] * sigma_ij12 * r6inv - lj4[itype][jtype] * sigma_ij6) -
          4.0 * epsilon[itype][jtype] * (ratio12 - ratio6);
  return factor_lj * philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJPoly::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0)
    return (void *)epsilon;
  return nullptr;
}
