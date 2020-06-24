# Naive Polycrystal Atomic Configuration Tool

Generate, naively, a polycrystalline sample using shortest Euclidean distance criteria ([i.e., Voronoi diagram](https://en.wikipedia.org/wiki/Voronoi_diagram)). The current version is capable of generating a LAMMPS configuration file given a specific crystal, number of grains, and box dimensions.

The [atomic simulation environment package (ASE)](https://wiki.fysik.dtu.dk/ase/index.html) is used for generating the crystallites. The user script just needs to supply the dictionary name from [materials.py](src/materials.py) to the '''Crystal''' class in [crystal.py](src/crystal.py). To add new materials/crystals you one just needs to provide a dictionary entry in [materials.py](src/materials.py) based on the other entry templates. The information to add new materials can be found in repositories such as the [materials project](https://materialsproject.org/) or [AFLOW](http://www.aflowlib.org/); in principal this could be done automatically by interfacing with those APIs but that would add another layer of dependencies.

This small code is essentially a rewrite of a [script](http://www.u.arizona.edu/~stefanb/Codes/graingen-v1.0.py) I used in graduate school. The idea is to generalize and marginally extend the capability.

## Dependencies

   - Numpy
   - Atomic Simulation Environment package

## Limitations

   - Periodic boundary conditions are not implemented so grain boundaries will always be present at box boundaries.
   - Only practical for small grains and simulation boxes.
   - Simulation box has to be orthorhombic.
