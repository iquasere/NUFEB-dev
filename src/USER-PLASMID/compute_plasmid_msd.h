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

#ifdef COMPUTE_CLASS

ComputeStyle(nufeb/plasmid/msd,ComputePlasmidMSD)

#else

#ifndef LMP_COMPUTE_PLASMID_MSD_H
#define LMP_COMPUTE_PLASMID_MSD_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePlasmidMSD : public Compute {
 public:
  ComputePlasmidMSD(class LAMMPS *, int, char **);
  virtual ~ComputePlasmidMSD();
  virtual void init();
  virtual void compute_vector();
  void set_arrays(int);

  class AtomVecBacillus *avec;

 protected:
  int nmsd;

  char *id_fix;
  class FixStore *fix_store;
  class FixPropertyPlasmid *fix_plasmid;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
