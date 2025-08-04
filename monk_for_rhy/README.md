#### Executable binaries

After compilation one can obtain the following binaries for photon propagation:

- `sphere`: spherical plasma cloud Comptonising primary radiation;
- `slab`: slab plasma cloud Comptonising primary radiation
- `ntdisc`: thermal photons from Novikov-Thorne disc;
- `3dcorona`: disc-extended corona system;

and the following three post-processing programs with which one can calculate
energy/polarisation spectrum or emissivity profile with files produced by the
programs above:

- `calspec`: calculate energy and polarisation spectrum
- `discen`: calculate energy of the reflecting photons as observed by disc fluid
- `emis`: calculate emissivity profile

For detailed descriptions, please run `doxygen doc.config` to generate the documentation file and look
for the links to the corresponding pages in `./html/files.html`. See below for a short introduction.

#### Building
##### System Requirements
- GNU Make
- c++ compiler that supports c++14 standard; for GCC starting from GCC 5.0
- openmp task feature; supported by GCC as of GCC 4.7
- libstdc++fs; included in GCC starting from GCC 5.3

##### Procedure
1. build `sim5` library under ./sim5 (I made some minor changes to sim5 source files, mostly header files, to be compatible with c++)
2. create `Makefile` from `Makefile.example` (I ignored `Makefile` in git as some variables are dependent on the user's system).
3. edit `Makefile` to assign the appropriate values to the following variables:
	- `SIM5INC`: path to directory that contains `sim5lib.h`
	- `SIM5OBJ`: path to `sim5lib.o`
	- `OBJDIR`: directory for object file
	- `BINDIR`: directory for binary files
	- `TRASHDIR`: directory for trash
	- `INSTALLDIR`: directory where you put executable binaries
4. edit `HOTXDIR` in `calhotcross.cpp`. This varaible should point to the directory that contains `logthetat.dat`, `logx.dat`, and `hotx.dat`. I put them under `./data`
5. `make objs` to build objects
6. `make binaries` to build the executable binaries
7. `make install` to cp the binaries to the desired directory for executable binaries (preferably a directory included in PATH0

#### Usage
##### Spherical/slab plasma cloud Comptonising primary radiation in flat space-time
1. Create a parameter file for `sphere/slab`
2. Run `sphere/slab`
3. Run `calspec` to calculate spectrum

##### Novikov-Thorne spectrum
1. Create a parameter file for `ntdisc`
2. Run `ntdisc`
3. Run `calspec` to calculate spectrum

##### AGN disc-corona system
1. Create a parameter file for `3dcorona` with `type = 0`;
2. Run `3dcorona` to produce `sca_params.dat` and `disk_params.dat`
3. Make a directory for disc photons; e.g., `disc`; and another directory for
   corona photons, e.g., `corona`
4. In `disc`, create an parameter file with `type = 1`, and run `3dcorona`
5. In `corona`, set the environment `OMP_NUM_THREADS`, create an parameter file with `type = 2`, and run `3dcorona`
6. In `disc` and `corona`, run `calspec` to calculate spectra of disc and corona, respectively
7. In `corona`, run `discen` to calculate energy of reflection photon in local frame, and then run `emissivity` to calculate the emissivity profile
