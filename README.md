Fortran METIS Interface
=======================

A modern Fortran interface to the METIS graph partitioning library

# Brief description

This is a Fortran interface to the [METIS software package](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) 
for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings
of sparse matrices. The interface makes use of the C interoperability features available in modern Fortran 
(i.e., Fortran 2003+) and provides a simple and safe way to call the original routines.

## METIS API

Example calls of the Fortran versions for the core METIS library routines:

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

## Object-oriented API

The object-oriented API is under development.

# Examples

Examples will be provided soon.

# Compiling

Compile instructions will be provided soon.

# Documentation

The latest API documentation can be found [here](). The documentation was generated from the source code using [FORD](https://github.com/cmacmackin/ford).

# License