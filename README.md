# tqchem
Source code of a class that builds a system of internal BAT coordinates for molecular systems.

# Installation
You can install the current version of tqchem via a conda install:
```
conda install tqchem -c conda-forge -c https://repo.prefix.dev/terraquantumag
```

# Running tqchem
After installation `tqchem` can be used as command line tool.
Two subcommands are provided:
1. `ttconf` to perform high accuracy conformer optimization
2. `generate` for a conventient way to
   convert smiles strings to .xyz files and generate some conformers quickly
   using RDkit.

To run TTconf an input structure needs to be provided either as smiles
string or as .xyz file
```sh
tqchem ttconf COCCNC
```
It is also possible to specify several starting conformers for the same molecule
as smiles string and .xyz files.
Further settings like the number of sweeps and the rank used in the
tensor train are summarized in the help message:
```sh
tqchem ttconf -h
```
