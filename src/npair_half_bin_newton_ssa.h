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

#ifdef NPAIR_CLASS

NPairStyle(half/bin/newton/ssa,
           NPairHalfBinNewtonSSA,
           NP_BIN | NP_NEWTON | NP_ORTHO | NP_SSA)

#else

#ifndef LMP_NPAIR_HALF_BIN_NEWTON_SSA_H
#define LMP_NPAIR_HALF_BIN_NEWTON_SSA_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairHalfBinNewtonSSA : public NPair {
 public:
  NPairHalfBinNewtonSSA(class LAMMPS *);
  ~NPairHalfBinNewtonSSA() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
