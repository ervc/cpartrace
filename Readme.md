# C-Partrace

- Author: Eric Van Clepper
- Version: 0.0.1

# Description
C implementation of python partrace model.

Main module in `partrace.c` and `partrace.h`, functions are kept in `src/` file. Make, by default, compiles and runs, can also call `make compile` or `make run` to do separately.

## TODO:

- [ ] Input files and reading in information
  - [x] Put alpha, aspect ratio, flaring angle, omegaframe, in the Model struct
  - [ ] Put the planet mass, location, and sun location in the model
  - [ ] read in particles from a file

- [ ] Grid stats while running
  - [ ] Residence times
  - [ ] velocities
    - [ ] put these in hdf5?

- [ ] Time dependence?
  - [ ] Make meshfields nt,nz,ny,nx?
  - [ ] Possibly include multiple planets

- [x] Integration
  - [x] Really fast!
  - [x] Full dt calculation
  - [x] Implement full coupling with gas before testing small particles
  - [x] Include diffusion
    - See bug about downward diffusion
  - [ ] Sometimes Fails after 8300 years. Sometime fails in the middle of writing to the file, which is weird...
    - This has not happened since the initial tests... Keep an eye out though

### Completed tasks
- Tests with rebound
  - Too complicated for now

- [x] Interpolation
  - [x]  when phi is between phi[nx-1] and phi[0]
  - [x]  When z is negative (when theta>pi/2)
  - [x]  when z is near the midplane (if `theta>theta[nz-1]` and `theta<pi-theta[nz-1]`)
  

## BUGS:

- [ ] Particle starts too fast, jumps outward in first year before slowing down
  - This is potentially not a bug, starting conditions are valid near midplane and I am starting them way up near the top of the disk...

### Fixed bugs
- [x] Particles mostly diffuse downward?
  - Needed to flip vz in the integration, veff was still using negative vz
- [x] Particles seem to spend more time in the negative z side than positive z
  - [x] I'm sure this has to do with where/when I'm checking for negative z values...
  -  I was not flipping the particle's velocity in determining the 
- [x] Theta out of range on particle number 6 in test. x0=y0=5au, z0=1.2au
  - this was not a bug. Theta really was out of the range of the simulation