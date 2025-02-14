.. NUFEB-manual documentation master file, created by
   sphinx-quickstart on Tue Oct 19 16:27:51 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NUFEB Documentation
========================================

NUFEB (Newcastle University Frontier in Engineering Biology) 
is an open source tool for 3D individual-based simulation of microbial communities. 
The tool is distributed under the terms of the GNU Public License (GPL). 

NUFEB is built on top of the molecular dynamic simulator LAMMPS (version: stable_29Oct2020), 
and extended with features for microbial modelling. LAMMPS is an open-source code, distributed under the terms of the
GPL. For more information about LAMMPS look at `LAMMPS homepage <https://lammps.sandia.gov/>`_.

----------

This site contains user and programmer documentation for the NUFEB use. If you find
errors or omissions in the manual or have suggestions for useful information to add, please
send an email to the developers.

| NUFEB Development team:
| Bowen Li, bowen.li2@newcastle.ac.uk
| Denis Taniguchi: denis.taniguchi@newcastle.ac.uk
| Joseph E. Weaver: joe.weaver@newcastle.ac.uk

----------


.. _user_documentation:
.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: User Guide
   :name: userdoc
   :includehidden:
   
   intro
   install
   install_win
   run_nufeb
   inputscript
   commands
   packages
   

.. _tutorial:
.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: Tutorial
   :name: progdoc
   :includehidden:

   work_microbe_soil
   tut_add_growth_model


.. _command_reference:
.. toctree::
   :name: reference
   :maxdepth: 1
   :caption: Command Reference

   list_general
   list_microbe_grid
   list_biology
   list_physics
   list_post_physics
   list_chemsitry
   list_post_chemsitry
   list_compute
   list_post
   list_run
   
   
******************
Indices and tables
******************

.. only:: html

   * :ref:`genindex`
   * :ref:`search`