# asevolution
## -- a cosmological and relativistic N-body code evolving the (a)symmetron
----------
## Compilation and usage

Before compilation, make sure that all required external libraries are
installed:

* LATfield2 (modified) (https://github.com/oyvach/LATfield2)
* FFTW version 3
* GNU Scientific Library (GSL) including CBLAS
* HDF5

Make sure that the include paths are set properly, or add them to the
makefile. Also check the compiler settings in the makefile. The code is
compiled by typing:

    make clean; make

A typical command to run a simulation looks like this:

    mpirun -np 16 ./asevolution -n 4 -m 4 -s settings.ini


For further information, please refer to the User Manual (manual.pdf)


Additional parameters added for AsGRD are explained in `metadata.hpp`


A settings.ini file that makes a cosmological run with the (a)symmetron may be found in ./settings/ together with a settings file for evolving the LCDM universe.

## animations
In order to use animation, use the LATfield2 version linked to above.
animation parameters are set in example ideal wall files.

## convention
* Snapshots and animations of the B-field are outputted with a factor $a^2\cdot$ sim.numpts too large and need to be adjusted in post-processing.
* Snapshots and animations of the hij_prime- and aq fields are outputted dimensionfully and need to be adjusted by division with e.g. H_prime in post-processing.

