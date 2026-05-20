# C-Partrace

- Author: Eric Van Clepper
- Version: 1.0.0

# Description
C implementation of python partrace model.

Main module in `partrace.c` and `partrace.h`, functions are kept in `src/` file. Use `make partrace_mpi` to make parallel implementation. To run, use `mpirun -np 4 ./partrace_mpi path/to/input.in`.

See [Van Clepper et al. 2025](https://iopscience.iop.org/article/10.3847/1538-4357/ada8a4) for model details.

## Fargo Integration

Partrace is originally designed to read in fargo outputs for the HD mesh.
Necessary files from FARGO are:

- domain_[x,y,z].dat
- variables.par
- planet0.dat*
- gas[dens, vx, vy, vz]N.dat**

*Current implementation assumes only 1 planet with data saved in
planet0.dat

**'N' is the output number from fargo. This is usually an int but can also be a
string, e.g., for averaged outputs.


## TODO:

- [ ] Input files and reading in information
  - [x] Put alpha, aspect ratio, flaring angle, omegaframe, in the Model struct
  - [x] Put the planet mass, location, and sun location in the model
  - [x] Read in from input file
  - [x] Read HD grid size from file
  - [x] Read in particles from a file

- [x] Parallelize main loop using open MPI
  - [ ] Merge outputs from different ranks into one file automatically
  - [ ] include separate build for non-parallel version?
  - [ ] Convert all `printf()` -> `if (rank == 0) printf()`
  - [ ] Fix error on build for `NLVL != 5`

- [ ] Time dependence
  - [ ] Make meshfields nt,nz,ny,nx?
  - [ ] Possibly include multiple planets

- [ ] Collisions
  - [ ] barnes hut tree
  - [ ] r* tree
  - [ ] relative velocities

### Completed tasks
- [x] Integration with other public HD models including
  - [x] FARGO3D ([Benítez-Llambay & Masset, 2016, ApJS, 223, 11](https://ui.adsabs.harvard.edu/abs/2016ApJS..223...11B/abstract))
  - [x] RADMC-like ([RADMC3D-2.0](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/))

- [x] Grid stats while running
  - [x] Residence times
  - [x] velocities

- [x] Make saving particle output optional
  - [x] Maybe make it so if fname=="NULL" then don't save

- [x] Interpolation
  - [x]  when phi is between phi[nx-1] and phi[0]
  - [x]  When z is negative (when theta>pi/2)
  - [x]  when z is near the midplane (if `theta>theta[nz-1]` and `theta<pi-theta[nz-1]`)
  

## BUGS:

Please contact me if you find any bugs in the current implementation!
