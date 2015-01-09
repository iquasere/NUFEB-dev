# -------- REQUIREMENTS: ---------
#  You must define your MOLTEMPLATE_PATH environment variable
#  and set it to the "common" subdirectory of your moltemplate distribution.
#  (See the "Installation" section in the moltemplate manual.)

# Create the coordinates of the atoms using PACKMOL
cd packmol_files

  packmol < mix_ethylene+benzene.inp
  mv -f system.xyz ../moltemplate_files/

cd ..



# Create LAMMPS input files this way:
cd moltemplate_files

  # Create the "oplsaa.lt" file which moltemplate will need

  cd oplsaa_lt_generator/
  ./oplsaa_moltemplate.py  oplsaa_subset.prm
  mv -f oplsaa.lt ..
  cd ..

  # run moltemplate

  moltemplate.sh -xyz system.xyz system.lt

  # This will generate various files with names ending in *.in* and *.data. 
  # Move them to the directory where you plan to run LAMMPS (in this case "../")
  mv -f system.data system.in* ../

  # Optional:
  # The "./output_ttree/" directory is full of temporary files generated by 
  # moltemplate.  They can be useful for debugging, but are usually thrown away.
  rm -rf output_ttree/

  # Optional:
  # Delete the "oplsaa.lt" file:
  rm -f oplsaa.lt


cd ../
