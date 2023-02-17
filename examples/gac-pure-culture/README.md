# Analysis of a methanogenic pure culture in presence of Granular Activated Carbon

## Run analysis

Run NUFEB simulations.
```
mpirun -np 1 ../../nufeb_mpi -in no_gac.lmp > no_gac.log
mpirun -np 1 ../../nufeb_mpi -in gac_sheet.lmp > sheet.log
mpirun -np 1 ../../nufeb_mpi -in gac_sinusoid.lmp > sinusoid.log
```

Parse output logs.
```
grep -E '^([^\s]+\s){9}[^\s]+$' NUFEB-dev/examples/gac_pure_culture/no_gac.log | tr -s ' ' '\t' | head -n -8  > no_gac.tsv
grep -E '^([^\s]+\s){9}[^\s]+$' NUFEB-dev/examples/gac_pure_culture/sinusoid.log | tr -s ' ' '\t' | head -n -8  > sinusoid.tsv
grep -E '^([^\s]+\s){9}[^\s]+$' NUFEB-dev/examples/gac_pure_culture/sheet.log | tr -s ' ' '\t' | head -n -8  > sheet.tsv
```

Python analysis. Requires pandas and matplotlib.
```
python gac_pure_culture.py
```

## Results

![Methane by GAC shape per number of cells](./examples/gac_pure_culture/plots/methane_by_number_cells.png)
