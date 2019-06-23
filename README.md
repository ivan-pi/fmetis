# Fortran METIS Interface

[![Build Status](https://travis-ci.com/ivan-pi/fmetis.svg?branch=master)](https://travis-ci.com/ivan-pi/fmetis)
[![Issues open](https://img.shields.io/github/issues/ivan-pi/fmetis.svg)](https://github.com/ivan-pi/fmetis/issues)

A modern Fortran interface to the [METIS software package](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings of sparse matrices.

* [Getting started](#getting-started)
    - [Specifying the METIS library location](#specifying-the-metis-library-location)
    - [Using a different compiler](#using-a-different-compiler)
    - [Specifying the integer and real type precision](#specifying-the-integer-and-real-type-precision)
    - [Choosing a custom install location](#choosing-a-custom-install-location)
* [Why Fortran METIS Interface?](#why-fortran-metis-interface)
* [METIS API](#metis-api)
* [Example usage](#example-usage)
* [Documentation](#documentation)
* [Contributing](#contributing)
* [License](#license)

## Getting started

To build and install METIS from source follow the instructions found [here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download). 

On Debian, Ubuntu and related Linux distributions you can install the METIS library, header file, and binaries through the terminal using:
```
apt-get install libmetis-dev libmetis5 metis
```
If using the [Homebrew](https://brew.sh/) package manager (e.g. on macOS or Linux) you can install METIS using:
```
brew install metis
```

Download the Fortran code:
```
git clone https://github.com/ivan-pi/fmetis
```

To build and test the Fortran METIS interface, type:
```
cd fmetis
mkdir build
cd build
cmake ..
make
ctest
```

Start using METIS in your code by including the following line in your program:
```Fortran
use metis_interface
```
To link your application with the module files and `fmetis` library in the `include/` and `lib/` folders (either in the local build folder or the installed location) use: `gfortran -o myapp -I path/to/include/ -L path/to/lib/ myapp.f90 -lfmetis -lmetis` 

Dependencies include:
* a Fortran 2018-compatible compiler,
* a working METIS library installation, and
* the CMake build system (v2.8 or higher).

### Specifying the METIS library location

For non-standard METIS install locations you can manually specify the right library location (either static or shared) like this:
```
cmake .. -DMETIS_LIB="custom/path/to/libmetis.a" # or libmetis.so
```

### Using a different compiler

If you want to build the Fortran interface with a different compiler (e.g. Intel Fortran, Flang, PGI) specify the `FC` environment variable when calling CMake:

```
FC=ifort cmake ..
```
Alternatively, you can set the CMake Fortran compiler flag the following way:
```
cmake .. -DCMAKE_Fortran_COMPILER=ifort
```
Note that the generated module files are not compatible between different compilers (or even different versions of a single compiler).

### Specifying the integer and real type precision

By default, the Fortran interface is built in single precision (32-bit integers and floating point number) matching the declarations in the original `metis.h` header file. Linking the 32 bit Fortran interface with a METIS library compiled using 64 bit integer precision will result in memory allocation errors at runtime. To configure the build for 64 bit signed integers, type:

```
cmake .. -DINT=64
```
A width of 64 should be specified if the number of vertices or the total
number of edges in the graph exceed the limits of a 32 bit signed integer
i.e., 2^31-1. Just like in the original METIS library, the `c_int32_t` and `c_int64_t` integer types used throughout the interface code must be supported by the Fortran compiler.

For double precision real variables, type:
```
cmake .. -DREAL=64
```

### Choosing a custom install location

To install the Fortran METIS interface in a non-default location (e.g. `~/.local/`, type:
```
cmake .. -DCMAKE_INSTALL_PREFIX="~/.local/"
make install
```

## Why Fortran METIS Interface?

Fortran is still one of the main programming languages used for scientific computing. In many finite element, finite volume, and other scientific codes, the fill-reducing orderings computed by METIS can help reduce the computational requirements of sparse matrix factorization by an order of magnitude. Moreover, the graph partitions produced by METIS can be used to divide meshes for parallel processing. 

While METIS was already designed with support for calls from Fortran, the C interoperability features available in modern Fortran (`>= 2003`) allow the METIS routines to be called in a simple and safe way that guarantess type checking and avoids issues with compiler name mangling.

This project was compiled with the aim to organize several resources on interfacing with METIS scattered throughout several internet forums, and thereby reduce the effort of calling METIS routines for new Fortran application developers.

## METIS API

The following core METIS library routines are available:

```Fortran
! Graph partitioning routines
ierr = METIS_ParthGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part)
ierr = METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part)

! Mesh partitioning routines
ierr = METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon,nparts,tpwgts,options,objval,epart,npart)
ierr = METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsize,nparts,tpwgts,options,objval,epart,npart)

! Sparse matrix reordering routines
ierr = METIS_NodeND(nvtxs,xadj,adjncy,vwgt,options,perm,iperm)

! Mesh-to-graph conversion routines
ierr = METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy)
ierr = METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy)

! Utility routines
ierr = METIS_SetDefaultOptions(options)
ierr = METIS_Free(ptr)
```

All functions return an integer status flag `ierr` that can be checked for success.

## Example usage

Examples of METIS usage from Fortran will be provided soon.

## Documentation

The latest documentation can be found [here](https://ivan-pi.github.io/fmetis/). The documentation was generated from the source code using [FORD](https://github.com/cmacmackin/ford). To generate the documentation locally run the following command in the root folder of the project (FORD must be installed first):
```
ford ./project.md
```
More information about METIS can be found on the [original homepage](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) of the [Karypis Lab](http://glaros.dtc.umn.edu/).
The original METIS documentation can be downloaded [here (PDF)](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).

## Contributing

Please submit a bug report or suggest an improvement by opening a [new issue](https://github.com/ivan-pi/fmetis/issues/new). Pull requests are also welcome and can be submitted [here](https://github.com/ivan-pi/fmetis/compare).

## License

The Fortran METIS Interface source code is being distributed under the [Apache License Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).  See the [LICENSE](https://raw.githubusercontent.com/ivan-pi/fmetis/master/LICENSE) file for further details. To succesfully use the Fortran METIS interface, the Fortran source code must be compiled and linked with the original METIS library that is distributed under the same Apache license (as of METIS version 5.0.3).

---

[Back to top](#fortran-metis-interface)
