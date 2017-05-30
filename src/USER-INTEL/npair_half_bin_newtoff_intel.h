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

NPairStyle(half/bin/newtoff/intel,
           NPairHalfBinNewtoffIntel,
           NP_HALF | NP_BIN | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_INTEL)

#else

#ifndef LMP_NPAIR_HALF_BIN_NEWTOFF_INTEL_H
#define LMP_NPAIR_HALF_BIN_NEWTOFF_INTEL_H

#include "npair_intel.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class NPairHalfBinNewtoffIntel : public NPairIntel {
 public:
  NPairHalfBinNewtoffIntel(class LAMMPS *);
  ~NPairHalfBinNewtoffIntel() {}
  void build(class NeighList *);

 private:
  template <class flt_t, class acc_t>
  void hbnni(NeighList *, IntelBuffers<flt_t,acc_t> *);
  template <class flt_t, class acc_t, int>
  void hbnni(const int, NeighList *, IntelBuffers<flt_t,acc_t> *, const int,
             const int);
};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
