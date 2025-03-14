CoordinationFortran
Overview

This repository contains different versions of a Fortran program for benchmarking, release, and debugging purposes.
Build Instructions

Use the following commands to build the respective versions:

    make - Build all versions.
    make benchmark - Build only the benchmarking version.
    make release - Build only the current release version.
    make debug - Build the debugging version (Ignore the long list of warnings, it's just code shaming me).

Configuration

Refer to setup.txt for most of the configuration details.
Notes

    There are currently two reading functions for lammps.data and trajectory.dmp, along with some redundant code that needs to be sorted and cleaned up.
    The lammps.data file is no longer technically needed as it was used in previous versions, but it should be sorted soon.
