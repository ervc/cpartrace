# C-Partrace

- Author: Eric Van Clepper
- Version: 1.0.0

# Description
C implementation of python partrace model.

Main module in `partrace.c` and `partrace.h`, functions are kept in `src/` file. Make, by default, compiles and runs, can also call `make compile` or `make run` to do separately.

See [Van Clepper et al. 2025](https://iopscience.iop.org/article/10.3847/1538-4357/ada8a4) for model details.

## TODO:

- [ ] Input files and reading in information
  - [x] Put alpha, aspect ratio, flaring angle, omegaframe, in the Model struct
  - [x] Put the planet mass, location, and sun location in the model
  - [x] Read in from input file
  - [ ] Read HD grid size from file
  - [ ] read in particles from a file

- [ ] Parallelize main loop using openMP

- [ ] Time dependence
  - [ ] Make meshfields nt,nz,ny,nx?
  - [ ] Possibly include multiple planets

- [ ] Collisions
  - [ ] barnes hut tree
  - [ ] r* tree

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

### Fixed bugs
- [x] Residence times is not correct because I'm gridding based on centers and not edges
- [x] Particles mostly diffuse downward?
  - Needed to flip vz in the integration, veff was still using negative vz
- [x] Particles seem to spend more time in the negative z side than positive z
  - [x] I'm sure this has to do with where/when I'm checking for negative z values...
  -  I was not flipping the particle's velocity in determining the 
- [x] Theta out of range on particle number 6 in test. x0=y0=5au, z0=1.2au
  - this was not a bug. Theta really was out of the range of the simulation
- [x] Particle starts too fast, jumps outward in first year before slowing down
  - Adjusted initial particle velocity to account for non-zero gas velocity