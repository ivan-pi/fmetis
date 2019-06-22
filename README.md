# Fortran METIS Interface

A modern Fortran interface to the [METIS software package](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings of sparse matrices.

* [Getting started](#getting-started)
* [Why Fortran METIS Interface?](#why-fortran-metis-interface)
* [METIS API](#metis-api)
* [Example usage](#example-usage)
* [Documentation](#documentation)
* [Contributing](#contributing)
* [License](#license)

## Getting started

To build and install METIS from source follow the instructions found [here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download). On Debian, Ubuntu and related Linux distributions you can install the METIS library, header file, and binaries using:
```
apt-get install libmetis-dev libmetis5 metis
```

Next, download and build the Fortran interface:
```
git clone https://github.com/ivan-pi/fmetis
cd fmetis
mkdir -p build
cd build
cmake ..
make
ctest
```

Use METIS in your code by including the following line:
```Fortran
use metis_interface
```

## Why Fortran METIS Interface?

Fortran is still one of the main programming languages used for scientific computing. In many finite element, finite volume, and other scientific codes, the fill-reducing orderings computed by METIS can help reduce the computational requirements of sparse matrix factorization by an order of magnitude. Moreover, the graph partitions produced by METIS can be used to divide meshes for parallel processing. 

While METIS was already designed with support for calls from Fortran, the C interoperability features available in modern Fortran (`>= 2003`) allow the METIS routines to be called in a simple and safe way with guaranteed type checking and avoiding issues with compiler name mangling.

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

Please submit a bug report or suggest a modification by opening a [new issue](https://github.com/ivan-pi/fmetis/issues/new). Pull requests are also welcome and can be submitted [here](https://github.com/ivan-pi/fmetis/compare).

## License

The Fortran METIS Interface source code is being distributed under the [Apache License Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).  See the [LICENSE](https://raw.githubusercontent.com/ivan-pi/fmetis/master/LICENSE) file for further details. To succesfully use the Fortran METIS interface, the Fortran source code must be [compiled and linked] with the original METIS library that is distributed under the same Apache license (as of METIS version 5.0.3). 

---

[Back to top](#fortran-metis-interface)