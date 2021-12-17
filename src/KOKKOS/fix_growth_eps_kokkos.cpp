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

#include <cstdio>
#include <cstring>
#include <cmath>
#include "fix_growth_eps_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "math_const.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixGrowthEPSKokkos<DeviceType>::FixGrowthEPSKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixGrowthEPS(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | MASK_MASK | RMASS_MASK | RADIUS_MASK | OUTER_MASS_MASK | OUTER_RADIUS_MASK;
  datamask_modify = RMASS_MASK | RADIUS_MASK | OUTER_MASS_MASK | OUTER_RADIUS_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixGrowthEPSKokkos<DeviceType>::update_cells()
{
  d_gmask = gridKK->k_mask.template view<DeviceType>();
  d_reac = gridKK->k_reac.template view<DeviceType>();
  d_dens = gridKK->k_dens.template view<DeviceType>();
  d_growth = gridKK->k_growth.template view<DeviceType>();

  gridKK->sync(execution_space, GMASK_MASK | CONC_MASK | DENS_MASK | REAC_MASK);

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixGrowthEPSCellsReactionTag>(0, grid->ncells), f);
  copymode = 0;

  gridKK->modified(execution_space, REAC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixGrowthEPSKokkos<DeviceType>::update_atoms()
{
  for (int i = 0; i < 3; i++) {
    boxlo[i] = domain->boxlo[i];
    grid_sublo[i] = grid->sublo[i];
    grid_subbox[i] = grid->subbox[i];
  }

  d_x = atomKK->k_x.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();
  d_outer_mass = atomKK->k_outer_mass.view<DeviceType>();
  d_outer_radius = atomKK->k_outer_radius.view<DeviceType>();

  d_conc = gridKK->k_conc.template view<DeviceType>();
  d_growth = gridKK->k_growth.template view<DeviceType>();
  d_gmask = gridKK->k_mask.template view<DeviceType>();

  cell_size = grid->cell_size;
  vol = cell_size * cell_size * cell_size;

  gridKK->sync(execution_space, GMASK_MASK | CONC_MASK | GROWTH_MASK);
  atomKK->sync(execution_space, datamask_read);

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixGrowthEPSCellsGrowthTag>(0, grid->ncells), f);

  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixGrowthEPSAtomsTag >(0, atomKK->nlocal), f);
  copymode = 0;

  atomKK->modified(execution_space, datamask_modify);
  gridKK->modified(execution_space, GROWTH_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixGrowthEPSKokkos<DeviceType>::Functor::Functor(FixGrowthEPSKokkos<DeviceType> *ptr):
  igroup(ptr->igroup), dt(ptr->dt), groupbit(ptr->groupbit),
  isub(ptr->isub), decay(ptr->decay),
  cell_size(ptr->cell_size), vol(ptr->vol),
  d_mask(ptr->d_mask), d_gmask(ptr->d_gmask), d_conc(ptr->d_conc), d_reac(ptr->d_reac),
  d_dens(ptr->d_dens), d_growth(ptr->d_growth), d_x(ptr->d_x), d_rmass(ptr->d_rmass),
  d_radius(ptr->d_radius), d_outer_mass(ptr->d_outer_mass), d_outer_radius(ptr->d_outer_radius)
{
  for (int i = 0; i < 3; i++) {
    boxlo[i] = ptr->boxlo[i];
    grid_sublo[i] = ptr->grid_sublo[i];
    grid_subbox[i] = ptr->grid_subbox[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixGrowthEPSKokkos<DeviceType>::Functor::operator()(FixGrowthEPSCellsReactionTag, int i) const
{
  if (!(d_gmask(i) & GHOST_MASK)) {
    d_reac(isub, i) += decay * d_dens(igroup ,i);
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixGrowthEPSKokkos<DeviceType>::Functor::operator()(FixGrowthEPSCellsGrowthTag, int i) const
{
  if (!(d_gmask(i) & GHOST_MASK)) {
    d_growth(igroup, i, 0) = -decay;
  }
}


/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixGrowthEPSKokkos<DeviceType>::Functor::operator()(FixGrowthEPSAtomsTag, int i) const
{
  if (d_mask[i] & groupbit) {
    // can't do:
    // int cell = grid->cell(d_x[i]);
    // hence copying the code here
    int c[3];
    const double small = 1e-12;

    for (int j = 0; j < 3; j++) {
      c[j] = static_cast<int>((d_x(i,j) - boxlo[j]) /
 			     cell_size + small) - grid_sublo[j];
    }

    int cell = c[0] + c[1] * grid_subbox[0] +
      c[2] * grid_subbox[0] * grid_subbox[1];

    static const double MY_PI  = 3.14159265358979323846; // pi
    const double density = d_rmass(i) /
      (4.0 * MY_PI / 3.0 * d_radius(i) * d_radius(i) * d_radius(i));

    d_rmass(i) = d_rmass(i) * (1 + d_growth(igroup, cell, 0) * dt);
    d_outer_mass(i) = 0.0;

    d_radius(i) = pow((3.0 / (4.0 * MY_PI)) * (d_rmass(i) / density), 1.0 / 3.0);
    d_outer_radius(i) = d_radius(i);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixGrowthEPSKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixGrowthEPSKokkos<LMPHostType>;
#endif
}
