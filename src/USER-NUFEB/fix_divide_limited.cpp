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

#include <iostream>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "lmptype.h"
#include "compute.h"
#include "math_const.h"
#include "random_park.h"
#include "modify.h"
#include "domain.h"
#include "atom_masks.h"

#include "lammps.h"
#include "region.h"
#include "fix_divide_limited.h"
#include "input.h"
#include "variable.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define DELTA 1.005

/* -------------------------Helper functions------------------------------ */

double putSphereOutsideFloorRegion(double sphere_z, double sphere_radius, double extent_zhi) {
  if (sphere_z - sphere_radius <= extent_zhi) {
    sphere_z = extent_zhi + sphere_radius;
  }
  return sphere_z;
}

void putSphereOutsideSinusoidalRegion(double (*sphere_coord)[3], double sphere_radius) {
  double surface_z = sin(sin(5e5 * (*sphere_coord)[0])) * sin(cos(5e5 * (*sphere_coord)[1])) / 10e4;
  if ((*sphere_coord)[2] - sphere_radius <= surface_z) {
    (*sphere_coord)[2] = surface_z + sphere_radius;
  }
}

/* ---------------------------------------------------------------------- */

FixDivideLimited::FixDivideLimited(LAMMPS *lmp, int narg, char **arg) :
    FixDivide(lmp, narg, arg)
{
  if (narg < 6)
    error->all(FLERR, "Illegal fix nufeb/division/limited command");

  diameter = utils::numeric(FLERR,arg[3],true,lmp);
  eps_density = utils::numeric(FLERR,arg[4],true,lmp);
  seed = utils::inumeric(FLERR,arg[5],true,lmp);

  int iarg = 6;
  nregion = -1;
  varflag = 0;
  vstr = xstr = ystr = zstr = nullptr;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region_blocked") == 0) {
      nregion = domain->find_region(arg[iarg+1]);
      if (nregion == -1) error->all(FLERR, "Fix nufeb/division/limited region ID does not exist");
      domain->regions[nregion]->init();
      domain->regions[nregion]->prematch();
      iarg += 2;
    } else if (strcmp(arg[iarg], "var") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nufeb/division/limited command");
      delete [] vstr;
      int n = strlen(arg[iarg+1]) + 1;
      vstr = new char[n];
      strcpy(vstr,arg[iarg+1]);
      varflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"set") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_atoms command");
      if (strcmp(arg[iarg+1],"x") == 0) {
        delete [] xstr;
        int n = strlen(arg[iarg+2]) + 1;
        xstr = new char[n];
        strcpy(xstr,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"y") == 0) {
        delete [] ystr;
        int n = strlen(arg[iarg+2]) + 1;
        ystr = new char[n];
        strcpy(ystr,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"z") == 0) {
        delete [] zstr;
        int n = strlen(arg[iarg+2]) + 1;
        zstr = new char[n];
        strcpy(zstr,arg[iarg+2]);
      } else error->all(FLERR,"Illegal create_atoms command");
      iarg += 3;
    } else {
      error->all(FLERR, "Illegal fix nufeb/division/limited command");
    }
  }

  // error check and further setup for variable test
  if (!vstr && (xstr || ystr || zstr))
    error->all(FLERR,"Incomplete use of variables in create_atoms command");
  if (vstr && (!xstr && !ystr && !zstr))
    error->all(FLERR,"Incomplete use of variables in create_atoms command");

  if (varflag) {
    vvar = input->variable->find(vstr);
    if (vvar < 0)
      error->all(FLERR,"Variable name for create_atoms does not exist");
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR,"Variable for create_atoms is invalid style");
    std::cout << "vvar: " << vvar << std::endl;
    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0)
        error->all(FLERR,"Variable name for create_atoms does not exist");
      if (!input->variable->internalstyle(xvar))
        error->all(FLERR,"Variable for create_atoms is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0)
        error->all(FLERR,"Variable name for create_atoms does not exist");
      if (!input->variable->internalstyle(yvar))
        error->all(FLERR,"Variable for create_atoms is invalid style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0)
        error->all(FLERR,"Variable name for create_atoms does not exist");
      if (!input->variable->internalstyle(zvar))
        error->all(FLERR,"Variable for create_atoms is invalid style");
    }
  }

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
}

/* ---------------------------------------------------------------------- */

FixDivideLimited::~FixDivideLimited()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void FixDivideLimited::compute()
{
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (atom->radius[i] * 2 >= diameter) {    // atom will divide
        double density = atom->rmass[i] / (     // mass
            4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]); // volume

        double split = 0.4 + (random->uniform() * 0.2);
        double imass = atom->rmass[i] * split;
        double jmass = atom->rmass[i] - imass;

        double iouter_mass = atom->outer_mass[i] * split;
        double jouter_mass = atom->outer_mass[i] - iouter_mass;

        double theta = random->uniform() * 2 * MY_PI;
        double phi = random->uniform() * (MY_PI);

        double oldx = atom->x[i][0];
        double oldy = atom->x[i][1];
        double oldz = atom->x[i][2];

        // update daughter cell i
        atom->rmass[i] = imass;
        atom->outer_mass[i] = iouter_mass;
        atom->radius[i] = pow(((6 * atom->rmass[i]) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
        atom->outer_radius[i] = pow((3.0 / (4.0 * MY_PI)) * ((atom->rmass[i] / density) + (iouter_mass / eps_density)), (1.0 / 3.0));

        double newx = oldx + (atom->radius[i] * cos(theta) * sin(phi) * DELTA);
        double newy = oldy + (atom->radius[i] * sin(theta) * sin(phi) * DELTA);
        double newz = oldz + (atom->radius[i] * cos(phi) * DELTA);

        double sphere_coord[3] = {newx, newy, newz};

        if (nregion != -1) {
          double floor_zhi = domain->regions[nregion]->extent_zhi;
          newz = putSphereOutsideFloorRegion(newz, atom->radius[i], floor_zhi);
        }
        if (varflag) {
          putSphereOutsideSinusoidalRegion(&sphere_coord, atom->radius[i]);
          newz = sphere_coord[2];
        }

        if (newx - atom->radius[i] < domain->boxlo[0]) {
          newx = domain->boxlo[0] + atom->radius[i];
        } else if (newx + atom->radius[i] > domain->boxhi[0]) {
          newx = domain->boxhi[0] - atom->radius[i];
        }
        if (newy - atom->radius[i] < domain->boxlo[1]) {
          newy = domain->boxlo[1] + atom->radius[i];
        } else if (newy + atom->radius[i] > domain->boxhi[1]) {
          newy = domain->boxhi[1] - atom->radius[i];
        }
        if (newz - atom->radius[i] < domain->boxlo[2]) {
          newz = domain->boxlo[2] + atom->radius[i];
        } else if (newz + atom->radius[i] > domain->boxhi[2]) {
          newz = domain->boxhi[2] - atom->radius[i];
        }
        atom->x[i][0] = newx;
        atom->x[i][1] = newy;
        atom->x[i][2] = newz;

        // create daughter cell j
        double jradius = pow(((6 * jmass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
        double jouter_radius = pow((3.0 / (4.0 * MY_PI)) * ((jmass / density) + (jouter_mass / eps_density)), (1.0 / 3.0));
        double *coord = new double[3];
        newx = oldx - (jradius * cos(theta) * sin(phi) * DELTA);
        newy = oldy - (jradius * sin(theta) * sin(phi) * DELTA);
        newz = oldz - (jradius * cos(phi) * DELTA);
        if (newx - jradius < domain->boxlo[0]) {
          newx = domain->boxlo[0] + jradius;
        } else if (newx + jradius > domain->boxhi[0]) {
          newx = domain->boxhi[0] - jradius;
        }
        if (newy - jradius < domain->boxlo[1]) {
          newy = domain->boxlo[1] + jouter_radius;
        } else if (newy + jradius > domain->boxhi[1]) {
          newy = domain->boxhi[1] - jradius;
        }
        if (newz - jradius < domain->boxlo[2]) {
          newz = domain->boxlo[2] + jradius;
        } else if (newz + jradius > domain->boxhi[2]) {
          newz = domain->boxhi[2] - jradius;
        }

        sphere_coord[0] = newx;
        sphere_coord[1] = newy;
        sphere_coord[2] = newz;

        if (nregion != -1) {
          double floor_zhi = domain->regions[nregion]->extent_zhi;
          newz = putSphereOutsideFloorRegion(newz, atom->radius[i], floor_zhi);
        }

        if (varflag) {
          putSphereOutsideSinusoidalRegion(&sphere_coord, atom->radius[i]);
          newz = sphere_coord[2];
        }

        coord[0] = newx;
        coord[1] = newy;
        coord[2] = newz;

        atom->avec->create_atom(atom->type[i], coord);
        int j = atom->nlocal - 1;

        atom->tag[j] = 0;
        atom->mask[j] = atom->mask[i];
        atom->v[j][0] = atom->v[i][0];
        atom->v[j][1] = atom->v[i][1];
        atom->v[j][2] = atom->v[i][2];
        atom->f[j][0] = atom->f[i][0];
        atom->f[j][1] = atom->f[i][1];
        atom->f[j][2] = atom->f[i][2];
        atom->omega[j][0] = atom->omega[i][0];
        atom->omega[j][1] = atom->omega[i][1];
        atom->omega[j][2] = atom->omega[i][2];
        atom->torque[j][0] = atom->torque[i][0];
        atom->torque[j][1] = atom->torque[i][1];
        atom->torque[j][2] = atom->torque[i][2];
        atom->rmass[j] = jmass;
        atom->biomass[j] = atom->biomass[i];
        atom->radius[j] = jradius;
        atom->outer_mass[j] = jouter_mass;
        atom->outer_radius[j] = jouter_radius;

        modify->create_attribute(j);

        for (int m = 0; m < modify->nfix; m++)
          modify->fix[m]->update_arrays(i, j);

        delete[] coord;
      }
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}
