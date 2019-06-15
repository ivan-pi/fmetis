! metis_interface.f90 -- Fortran METIS Interface
!
! Copyright (C) 2018 Ivan Pribec <ivan.pribec@gmail.com>
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.


!*****************************************************************************************
!> author: Ivan Pribec
!  date: 7/2018
!  license: MIT
!
! A Fortran interface to the METIS graph partitioning library.
!
module metis_interface

    use iso_c_binding, only : c_int, c_double, c_ptr
    
    implicit none
    private

    ! Graph partitioning routines
    public :: METIS_PartGraphRecursive
    public :: METIS_PartGraphKway

    ! Mesh partitioning routines
    public :: METIS_PartMeshDual
    public :: METIS_PartMeshNodal

    ! Sparse matrix reordering routines
    public :: METIS_NodeND

    ! Mesh-to-graph conversion routines
    public :: METIS_MeshToDual
    public :: METIS_MeshToNodal

    ! Utility routines
    public :: METIS_SetDefaultOptions
    public :: METIS_Free

    ! Constants
    integer, parameter, public :: METIS_NOPTIONS = 40 !! Number of METIS options.


    !
    ! METIS' API
    !

    ! http://glaros.dtc.umn.edu/gkhome/node/877

    interface

    !
    ! Graph partitioning routines
    !

!*****************************************************************************************
!> This function is used to partition a graph into `nparts` parts using recursive bisection.
! 
!  If `tpwgt` is present, the *target partition weight* for the `i`-th partition and `j`-th constraint should
!  be specified at `tpwgts(i*ncon+j)` (the numbering for both partitions and constraints starts from 0).
!  For each constraint, the sum of the `tpwgts`entries must be 1.0.
!
!  The following options are valid: <br />
!  `METIS_OPTION_CTYPE`, `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, 
!  `METIS_OPTION_NO2HOP`, `METIS_OPTION_NCUTS`, `METIS_OPTION_NITER`, 
!  `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`, `METIS_OPTION_NUMBERING`,
!  `METIS_OPTION_DBGLVL`
!
    function METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,&
        vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part) result(ierr) bind(C,name="METIS_PartGraphRecursive")
        import c_int, c_double, METIS_NOPTIONS
        
        ! Parameters
        integer(c_int), intent(in) :: nvtxs
            !! The number of vertices in the graph.
        integer(c_int), intent(in) :: ncon
            !! The number of balancing constraints. It should be atleast 1.
        integer(c_int), intent(in), dimension(*) :: xadj, adjncy
            !! The adjacency structure of the graph as described in section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), dimension(*), optional :: vwgt ! NULL
            !! The weights of the vertices as described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), dimension(*), optional :: vsize ! NULL
            !! The size of the vertices for computing the total communication volume as described in section 5.7 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), dimension(*), optional :: adjwgt ! NULL
            !! The weights of the edges as describe in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int) , intent(in) :: nparts
            !! The number of parts to partition the graph.
        real(c_double), intent(in), optional :: tpwgts(nparts*ncon)
            !! An array of size `nparts*ncon` that specifies the desired weight for each partition and constraint.
            !! If not present, the graph is divided equally among the partitions. More in the description.
        real(c_double), intent(in), optional :: ubvec(ncon)
            !! An array of size `ncon` that specifies the allowed load imbalance for each constraint. 
            !! For the `i`-th partition and `j`-th constraint the allowed weight is the `ubvec(j)*tpwgts(i*ncon+j)`
            !! fraction of the `j`-th's constraint total weight. If not present, the load imbalance
            !! tolerance is 1.001 (for `ncon = 1`) or 1.01 (for `ncon > 1`).
        integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the METIS manual. See description for valid options.
        integer(c_int), intent(out) :: objval
            !! Upon successful completion, this variable stores the edge-cut or the total communication volume of the partitioning
            !! solution. The value returned depends on the partitioning's objective function.
        integer(c_int), intent(out) :: part(nvtxs)
            !! This is a vector of size `nvtxs` that upon successful completion stores the partition vector of the graph.
            !! The numbering of this vector starts from either 0 or 1, depending on the value of `options(METIS_OPTION_NUMBERING)`.

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> This function is used to partition a graph into `nparts` parts using multilevel \(k\)-way partitioning.
!
!  If `tpwgt` is present, the *target partition weight* for the `i`-th partition and `j`-th constraint should
!  be specified at `tpwgts(i*ncon+j)` (the numbering for both partitions and constraints starts from 0).
!  For each constraint, the sum of the `tpwgts`entries must be 1.0.
!
!  The following options are valid: <br />
!  `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`, `METIS_OPTION_IPTYPE`,
!  `METIS_OPTION_RTYPE`, `METIS_OPTION_NO2HOP`, `METIS_OPTION_NCUTS`,
!  `METIS_OPTION_NITER`, `METIS_OPTION_UFACTOR`, `METIS_OPTION_MINCONN`,
!  `METIS_OPTION_CONTIG`, `METIS_OPTION_SEED`, `METIS_OPTION_NUMBERING`,
!  `METIS_OPTION_DBGLVL`
!
    function METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,&
        vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part) result(ierr) bind(C,name="METIS_PartGraphKway")
        import c_int, c_double, METIS_NOPTIONS
        
        ! Parameters
        integer(c_int), intent(in) :: nvtxs
            !! The number of vertices in the graph.
        integer(c_int), intent(in) :: ncon
            !! The number of balancing constraints. It should be at least 1.
        integer(c_int), intent(in), dimension(*) :: xadj, adjncy
            !! The adjacency structure of the graph as described in section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), dimension(*), optional :: vwgt
            !! The weights of the vertices as described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), dimension(*), optional :: vsize
            !! The size of the vertices for computing the total communication volume as described in section 5.7 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), dimension(*), optional :: adjwgt
            !! The weights of the edges as describe in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in) :: nparts
            !! The number of parts to partition the graph.
        real(c_double), intent(in), optional :: tpwgts(nparts*ncon)
            !! An array of size `nparts*ncon` that specifies the desired weight for each partition and constraint.
            !! If not present, the graph is divided equally among the partitions. More in the description.
        real(c_double), intent(in), optional :: ubvec(ncon)
            !! An array of size `ncon` that specifiew the allowed load imbalance for each constraint. 
            !! For the `i`-th partition and `j`-th constraint the allowed weight is the `ubvec(j)*tpwgts(i*ncon+j)`
            !! fraction of the `j`-th's constraint total weight. If not present, the load imbalance
            !! tolerance is 1.001 (for `ncon == 1`) or 1.01 (for `ncon > 1`).
        integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(c_int), intent(out) :: objval
            !! Upon successful completion, this variable stores the edge-cut or the total communication volume of the partitioning
            !! solution. The value returned depends on the partitioning's objective function.
        integer(c_int), intent(out) :: part(nvtxs)
            !! This is a vector of size `nvtxs` that upon successful completion stores the partition vector of the graph.
            !! The numbering of this vector starts from either 0 or 1, depending on the value of `options(METIS_OPTION_NUMBERING)`.
        
        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Mesh partitioning routines
    !

!*****************************************************************************************
!> This function is used to partition a mesh into `nparts` parts based on a partitioning of the mesh's dual graph.
!
!  The following options are valid: <br />
!  `METIS_OPTION_PTYPE`, `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`,
!  `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NCUTS`,
!  `METIS_OPTION_NITER`, `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`,
!  `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
!
    function METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon, &
        nparts,tpwgts,options,objval,epart,npart) result(ierr) bind(C,name="METIS_PartMeshDual")
        import c_int, METIS_NOPTIONS
        integer(c_int), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(c_int), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(c_int), intent(in), dimension(*) :: eptr,eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), optional :: vwgt(ne)
            !! An array of size `ne` specifying the weights of the elements. If not present,
            !! all elements have an equal weight.
        integer(c_int), intent(in), optional :: vsize(ne)
            !! An array of size `ne` specifying the size of the elements that is used
            !! for computing the total comunication volume as described in Section 5.7 of the  [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
            !! If not present, the objective is cut or all elements have an equal size.
        integer(c_int), intent(in) :: ncommon
        integer(c_int), intent(in) :: nparts
            !! The number of parts to partition the mesh.
        integer(c_int), intent(in), optional :: tpwgts(nparts)
            !! An array of size `nparts` that specifies the desired weight for each partition. The *target
            !! partition weight* for the `i`-th partition is specified at `tpwgts(i)` (the numbering for the 
            !! partitions starts from 0). The sum of the `tpwgts` entries must be 1.0. <br /> If not present, the graph
            !! is divided equally among the partitions.
        integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(c_int), intent(out) :: objval
            !! Upon successful completion, this variable stores either the edgecut or the total communication
            !! volume of the dual graph's partitioning.
        integer(c_int), intent(out) :: epart(ne)
            !! A vector of size `ne` that upon successful completion stores the partition vector for the elements
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.
        integer(c_int), intent(out) :: npart(nn)
            !! A vector of size `nn` that upon successful completion stores the partition vector for the nodes 
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> This function us used to partition a mesh into `nparts` parts based on a 
!  partitioning of the mesh's nodal graph.
!
!  The following options are valid: <br />
!  `METIS_OPTION_PTYPE`, `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`,
!  `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NCUTS`,
!  `METIS_OPTION_NITER`, `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`,
!  `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
!
    function METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsize, &
        nparts,tpwgts,options,objval,epart,npart) result(ierr) bind(C,name="METIS_PartMeshNodal")
        import c_int, METIS_NOPTIONS
        integer(c_int), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(c_int), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(c_int), intent(in), dimension(*) :: eptr,eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), optional :: vwgt(nn)
            !! An array of size `nn` specifying weights of the nodes. If not passed, all nodes have an equal weight.
        integer(c_int), intent(in), optional :: vsize(nn)
            !! An array of size `nn` specifying the size of the nodes that is used for computing the
            !! total comunication volume as described in Section 5.7 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). If not passed,
            !! the objective is cut or all nodes have an equal size.
        integer(c_int), intent(in) :: nparts
            !! The number of parts to partition the mesh.
        integer(c_int), intent(in), optional :: tpwgts(nparts)
            !! An array of size `nparts` that specifies the desired weight for each partition. The *target
            !! partition weight* for the `i`-th partition is specified at `tpwgts(i)` (the numbering for the 
            !! partitions starts from 0). The sum of the `tpwgts` entries must be 1.0. If not passed, the graph
            !! is divided equally among the partitions.
        integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(c_int), intent(out) :: objval
            !! Upon successful completion, this variable stores either the edgecut or the total communication
            !! volume of the nodal graph's partitioning.
        integer(c_int), intent(out) :: epart(ne)
            !! A vector of size `ne` that upon successful completion stores the partition vector for the elements
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.
        integer(c_int), intent(out) :: npart(nn)
            !! A vector of size `nn` that upon successful completion stores the partition vector for the nodes 
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Sparse Matrix Reordering Routines
    !

!*****************************************************************************************
!> This function computes fill reducing orderings of sparse matrices using the
!  multilevel nested dissection algorithm.
!
!  Let \(A\) be the original matrix and \(A'\) be the permuted matrix. The arrays
!  `perm` and `iperm` are defined as follows. Row (column) `i` of \(A'\) is the
!  `perm(i)` row (column) of \(A\), and row (column) `i` of \(A\) is the `iperm(i)`
!  row (column) of \(A'\). the numbering of this vector starts from either 0 or 1,
!  depending on the value of `options(METIS_OPTION_NUMBERING)`.
!
!  If the graph is weighted, meaning `vgwt` was provided, the nested dissection ordering computes
!  vertex separators that minimize the sum of the weights of the vertices on the separators.
!
!  The following options are valid: <br />
!  `METIS_OPTION_CTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NO2HOP`,
!  `METIS_OPTION_NSEPS`, `METIS_OPTION_NITER`, `METIS_OPTION_UFACTOR`,
!  `METIS_OPTION_COMPRESS`, `METIS_OPTION_CCORDER`, `METIS_OPTION_SEED`,
!  `METIS_OPTION_PFACTOR`, `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
!
!#Example
! The code below generates the nodal graph of the following mesh:
!```Fortran
!integer(c_int), parameter :: n = 15, m = 22
!integer(c_int) :: xadj(n+1), adjncy(2*m)
!integer(c_int) :: perm(n), iperm(n)
!integer(c_int) :: options(0:39), ierr
!
!xadj = [1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45]
!adjncy = [2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6, &
!                  8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,7,11,13, &
!                  8,12,14,9,13,15,10,14]
!
!ierr = METIS_SetDefaultOptions(options)
!options(18) = 1
!
!ierr = METIS_NodeND(n,xadj,adjncy,options=options,perm=perm,iperm=iperm)
!end
!```
    function METIS_NodeND(nvtxs,xadj,adjncy,vwgt,options,perm,iperm) result(ierr) bind(C,name="METIS_NodeND")
        import c_int, METIS_NOPTIONS

        ! Parameters
        integer(c_int), intent(in) :: nvtxs
            !! The number of vertices in the graph.
        integer(c_int), intent(in), dimension(*) :: xadj, adjncy
            !! The adjacency structure of the graph as described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in), optional :: vwgt(nvtxs)
            !! An array of size `nvtxs` specifying the weights of the vertices.
        integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            !! This is the array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(c_int), intent(out) :: perm(nvtxs), iperm(nvtxs)
            !! Vectors of size `nvtxs`. Upon successful completion, they store the fill-reducing
            !! permutation and inverse-permutation. More in the description.

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Mesh-to-graph conversion routines
    !

!*****************************************************************************************
!> This function is used to generate the dual graph of a mesh.
!
!@note
! To use the returned arrays `xadj` and `adjncy`, these must be first converted from
! a C pointer to a Fortran pointer using the subroutine `c_f_pointer(cptr,fptr,shape)` 
! that assigns the target of the C pointer `cptr` to the Fortran pointer `fptr` and
! specifies its shape. The `shape` is an integer rank-one array, storing the size `ne+1`
! in case of the dual graph. The size of the new `adjncy` array is stored in the 
! last element of `xadj` when using C-style numbering. An example is shown below.
!@endnote
!
!@warning
! Memory for the returned arrays `xadj` and `adjncy` is allocated by METIS' API in C
! using the standard `malloc` function. It is the responsibility of the application to free
! this memory by calling `free`. Therefore, METIS provides the [[METIS_Free]] function that is a wrapper to
! C's `free`function.
!@endwarning
!
!# Example
! The code below generates the nodal graph of the following mesh:
! __image__
!```Fortran
!use iso_c_binding, only : c_int, c_ptr, c_f_pointer
!integer(c_int), parameter :: ne = 3, nn = 8, npel = 4
!integer(c_int) :: ierr, eptr(ne+1), eind(ne*npel), numflag, ncommon
!type(c_ptr) :: xadj, adjncy
!integer(c_int), dimension(:), pointer :: fxadj => null(), fadjncy => null()
!
!numflag = 0 ! C-style numbering
!ncommon = 2 ! 2 common nodes per edge
!eptr = [0,4,8,12]
!eind = [0,1,2,3,1,4,5,2,4,6,7,5] ! note four nodes per element
!
!ierr = METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy)
!call c_f_pointer(xadj,fxadj,shape=[ne+1]) ! xadj is of the size ne+1
!call c_f_pointer(adjncy,fadjncy,shape=[fxadj(ne+1)]) ! last value in xadj is the size of adjncy
!
! !... use values in fxadj and fadjncy ...
!
!call METIS_Free(xadj)
!call METIS_Free(adjncy)
!end
!```
    function METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy) result(ierr) bind(C,name="METIS_MeshToDual")
        import c_int, c_ptr
        
        ! Parameters
        integer(c_int), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(c_int), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(c_int), intent(in), dimension(*) :: eptr, eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in) :: ncommon
            !! The number of common nodes that two elements must have in order to put
            !! an edge between them in the dual graph.
        integer(c_int), intent(in) :: numflag
            !! Used to indicate which numbering scheme is used for `eptr` and `eind`. 
            !! The possible values are: <br />
            !! 0 - C-style numbering is assumed that starts from 0 <br />
            !! 1 - Fortran-style numbering is assumed that starts from 1
        type(c_ptr), intent(out) :: xadj, adjncy
            !! These arrays store the adjacency structure of the generated dual graph. 
            !! The format of the adjacency structure is described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> This function is used to generate the nodal graph of a mesh.
!
!@note
! To use the returned arrays `xadj` and `adjncy`, these must be first converted from
! a C pointer to a Fortran pointer using the subroutine `c_f_pointer(cptr,fptr,shape)` 
! that assigns the target of the C pointer `cptr` to the Fortran pointer `fptr` and
! specifies its shape. The `shape` is an integer rank-one array, storing the size `nn+1`
! in case of the nodal graph. The size of the new `adjncy` array is stored in the 
! last element of `xadj` when using C-style numbering. An example is shown below.
!@endnote
!
!@warning
! Memory for the returned arrays `xadj` and `adjncy` is allocated by METIS' API in C
! using the standard `malloc` function. It is the responsibility of the application to free
! this memory by calling `free`. Therefore, METIS provides the [[METIS_Free]] function that is a wrapper to
! C's `free`function.
!@endwarning
!
!# Example
! The code below generates the nodal graph of the following mesh:
! __image__
!```Fortran
!use iso_c_binding, only : c_int, c_ptr, c_f_pointer
!integer(c_int), parameter :: ne = 3, nn = 8, npel = 4
!integer(c_int) :: ierr, eptr(ne+1), eind(ne*npel), numflag
!type(c_ptr) :: xadj, adjncy
!integer(c_int), dimension(:), pointer :: fxadj => null(), fadjncy => null()
!
!numflag = 0 ! C-style numbering
!eptr = [0,4,8,12]
!eind = [0,1,2,3,1,4,5,2,4,6,7,5] ! note four nodes per element
!
!ierr = METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy)
!call c_f_pointer(xadj,fxadj,shape=[nn+1]) ! xadj is of the size nn+1
!call c_f_pointer(adjncy,fadjncy,shape=[fxadj(nn+1)]) ! last value in xadj is the size of adjncy
!
! !... use values in fxadj and fadjncy ...
!
!call METIS_Free(xadj)
!call METIS_Free(adjncy)
!end
!```
    function METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy) result(ierr) bind(C,name="METIS_MeshToNodal")
        import c_int, c_ptr
        
        ! Parameters
        integer(c_int), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(c_int), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(c_int), intent(in), dimension(*) :: eptr, eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(c_int), intent(in) :: numflag
            !! Used to indicate which numbering scheme is used for `eptr` and `eind`.
            !! The possible values are: <br />
            !! 0 - C-style numbering is assumed that starts from 0 <br />
            !! 1 - Fortran-style numbering is assumed that starts from 1
        type(c_ptr), intent(out) :: xadj, adjncy
            !! These arrays store the adjacency structure of the generated dual graph. 
            !! The format of the adjacency structure is described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Utility routines
    !

!*****************************************************************************************
!> Initializes the options array into its default values.
!
!@note
! The passed array `options` must have the size `METIS_NOPTIONS` (40).
! To be able to use the parameters in the [[metis_enum]] module it is recommended to use
! zero-based indexing for the options array: 
!```Fortran
!integer(c_int) :: opts(0:39)
!```
!@endnote
!
!# Examples
! To set Fortran style index-numbering use:
!```Fortran
!integer :: opts(0:39)
! 
!call METIS_SetDefaultOptions(opts)
!opts(17) = 1 ! Fortran-style index numbering
!```
!
! Other options can be changed using parameters from the [[metis_enum]] module.
!```Fortran
!use metis_interface, only : METIS_SetDefaultOptions
!use metis_enum, only : METIS_OPTION_DBGLVL, METIS_DBG_INFO
!integer :: opts(0:39)
! 
!call METIS_SetDefaultOptions(opts)
!opts(METIS_OPTION_DBGLVL) = METIS_DBG_INFO ! Show various diagnostic messages
!end
!```
    function METIS_SetDefaultOptions(options) result(ierr) bind(C,name="METIS_SetDefaultOptions")
        import c_int, METIS_NOPTIONS
        
        ! Parameters
        integer(c_int), intent(out) :: options(METIS_NOPTIONS)
            !! The array of options that will be initialized.

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> Frees the memory that was allocated by either the [[METIS_MeshToDual]] or the
!  [[METIS_MeshToNodal]] routines for returning the dual or nodal graph of a mesh.
!
!@warning Memory deallocation should always happen on the same side it was allocated!
! Also check the descriptions of the above-mentioned routines.
!
!# Example
!
!```Fortran
! type(c_ptr) :: xadj(:),adjncy(:)
! 
! call METIS_MeshToNodal(...,xadj,adjncy)
! 
! ! xadj and adjncy should be deallocated on the C side! ;)
! call METIS_Free(xadj)
! call METIS_Free(adjncy)
!```   
    function METIS_Free(ptr) result(ierr) bind(C,name="METIS_Free")
        import c_int, c_ptr
        
        ! Parameters
        type(c_ptr), value :: ptr
            !! The pointer to be freed. This pointer should be one of the `xadj` or `adjncy`
            !! arrays returned by METIS' API routines.

        ! Returns
        integer(c_int) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.
    end function
!*****************************************************************************************

    end interface

end module
!*****************************************************************************************





