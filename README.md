# Fortran METIS Interface

A modern Fortran interface to the [METIS software package](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings of sparse matrices.

* [Getting started](#getting-started)
* [Why Fortran METIS Interface?](#why-fortran-metis-interface)
* [METIS API](#metis-api)
* [Object-oriented API](#object-oriented-api)
* [Example usage](#example-usage)
* [Documentation](#dcoumentation)
* [Contributing](#contributing)
* [License](#license)

## Getting started

Download and install METIS according to instructions found [here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download).

Download and build the Fortran interface:

```
git clone https://github.com/ivan-pi/fmetis
cd fmetis
mkdir -p build
cd build
cmake ..
make
ctest
```

Use METIS in your code by including the modules:

```Fortran
use metis_enum
use metis_interface
```

## Why Fortran METIS Interface?

Fortran is still one of the main programming languages used for scientific computing. In many finite element, finite volume, and other scientific codes, the fill-reducing orderings computed by METIS can help reduce the computational requirements of sparse matrix factorization. Moreover, the graph partitions produced by METIS can be used to divide meshes for parallel processing. 

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

## Object-oriented API

The object-oriented API is still under development.

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

Please submit a bug report or suggest a modification by opening a [new issue](https://github.com/ivan-pi/fmetis/issues/new).

## License

The Fortran METIS Interface source code, related files and documentation are distributed under the permissive MIT software license.  See the [LICENSE](https://raw.githubusercontent.com/ivan-pi/fmetis/master/LICENSE) file for more details. Note that to succesfully use the Fortran METIS interface, the source code of this project must be [compiled and linked](#compiling-and-linking) with the original METIS source files that are distributed under the [Apache License Version 2.0](http://www.apache.org/licenses/LICENSE-2.0). 

## To-Do List
- [ ] Provide usage examples for partitioning graphs, partitioning meshes and computing sparse matrix reorderings
- [ ] Provide examples for reading and writing METIS graph and mesh files
- [ ] Move documentation to a separate branch