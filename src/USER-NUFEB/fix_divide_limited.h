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

#ifdef FIX_CLASS

FixStyle(nufeb/division/limited,FixDivideLimited)

#else

#ifndef LMP_FIX_DIVIDE_LIMITED_H
#define LMP_FIX_DIVIDE_LIMITED_H

#include "fix_divide.h"

double putSphereOutsideFloorRegion(double sphere_z, double sphere_radius, double extent_zhi);
void replaceVariables(std::string* expression, double x_val, double y_val);
void putSphereOutsideSinusoidalRegion(double (*sphere_coord)[3], double sphere_radius);

namespace LAMMPS_NS {

  class FixDivideLimited : public FixDivide {
  public:
    FixDivideLimited(class LAMMPS *, int, char **);
    virtual ~FixDivideLimited();
    virtual void compute();

  private:
    int varflag,vvar,xvar,yvar,zvar;
    char *vstr,*xstr,*ystr,*zstr;

  protected:
    double diameter;
    double eps_density;
    int seed;

    int nregion;
    int region_blocked;

    class RanPark *random;
  };

}

#endif
#endif

/* ERROR/WARNING messages:
*/
