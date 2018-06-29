set(STANDARD_PACKAGES ASPHERE BODY CLASS2 COLLOID COMPRESS CORESHELL DIPOLE GPU
                      GRANULAR KIM KOKKOS KSPACE LATTE MANYBODY MC MEAM MISC
                      MOLECULE MPIIO MSCG OPT PERI POEMS
                      PYTHON QEQ REAX REPLICA RIGID SHOCK SNAP SRD VORONOI)

set(USER_PACKAGES USER-ATC USER-AWPMD USER-BOCS USER-CGDNA USER-CGSDK USER-COLVARS
                  USER-DIFFRACTION USER-DPD USER-DRUDE USER-EFF USER-FEP USER-H5MD
                  USER-INTEL USER-LB USER-MANIFOLD USER-MEAMC USER-MESO
                  USER-MGPT USER-MISC USER-MOFFF USER-MOLFILE
                  USER-NETCDF USER-OMP USER-PHONON USER-QMMM USER-QTB
                  USER-QUIP USER-REAXC USER-SMD USER-SMTBQ USER-SPH USER-TALLY
                  USER-UEF USER-VTK)

set(PACKAGES_WITH_LIB COMPRESS GPU KIM KOKKOS LATTE MEAM MPIIO MSCG POEMS PYTHON REAX VORONOI
                      USER-ATC USER-AWPMD USER-COLVARS USER-H5MD USER-LB USER-MOLFILE
                      USER-NETCDF USER-QMMM USER-QUIP USER-SMD USER-VTK)

set(ALL_PACKAGES ${STANDARD_PACKAGES} ${USER_PACKAGES})

foreach(PKG ${USER_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
