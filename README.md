# Fortran METIS Interface

## Brief description

This package provides a Fortran interface to the [METIS software package](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings of sparse matrices. The interface makes use of the C interoperability features available in modern Fortran (i.e., Fortran 2003+) and provides a simple and safe way to call the original routines.

### METIS API

Examples of calls to the Fortran versions of the core METIS library routines:

```Fortran
! All functions return an integer status flag ierr

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

### Object-oriented API

The object-oriented API is currently under development.

## Examples

Examples of METIS usage from Fortran will be provided soon.

## Compiling and linking

Compile instructions will be provided soon.

## Documentation

The latest documentation can be found [here](https://ivan-pi.github.io/fmetis/). The documentation was generated from the source code using [FORD](https://github.com/cmacmackin/ford). To generate the documentation locally run the following command in the root folder of the project (FORD must be installed first):
```
ford ./project.md
```
More information about METIS can be found on the [original homepage](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) of the [Karypis Lab](http://glaros.dtc.umn.edu/).
The original METIS documentation may be downloaded [here (PDF)](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).

## License

The Fortran METIS Interface source code, related files and documentation are distributed under the permissive MIT software license.  See the [LICENSE](https://raw.githubusercontent.com/ivan-pi/fmetis/master/LICENSE) file for more details. Note that to succesfully use the Fortran METIS interface, the source code of this project must be [compiled and linked](#compiling-and-linking) with the original METIS source files that are distributed under the [Apache License Version 2.0](http://www.apache.org/licenses/LICENSE-2.0). 

## To-Do List
- [ ] Provide usage examples for partitioning graphs, partitioning meshes and computing sparse matrix reorderings
- [ ] Provide examples for reading and writing METIS graph and mesh files
- [ ] Move documentation to a separate branch
- [ ] Develop unit tests for the "raw" C and object-oriented interfaces
