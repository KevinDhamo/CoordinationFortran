make to get all versions
make benchmark - for just the benchmarking version
make release   - for just the current version
make debug     - for debugging version
Ignore the long list of warnings from debug version, its just code shaming me.

setup.txt contains most of configuration things
Currently there are two reading functions for lammps.data and for trajectory.dmp along with redundant code about them that need 
to be sorted and cleaned up. 
It doesn't technically need the lammps.data file anymore as it was used in previous version but it should be sorted soon.
