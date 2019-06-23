!>  Example of partitioning a mesh composed of 4 quadrilaterals and 9 nodes
!   based upon its nodal graph.
!
!   Source: https://www.cfd-online.com/Forums/main/112366-using-metis-functions-fortran.html#post404734
!
program test_PartMeshNodal2

    use metis_interface, only: idx_t, METIS_SetDefaultOptions, METIS_PartMeshNodal, METIS_NOPTIONS, &
        METIS_OPTION_NUMBERING, METIS_OPTION_CONTIG, METIS_OK
    implicit none

    integer(idx_t), parameter :: ne = 4    ! number of elements
    integer(idx_t), parameter :: nn = 9    ! number of nodes

    integer(idx_t) :: eptr(ne+1), eind(4*ne)   ! arrays storing mesh structure
    integer(idx_t) :: epart(ne), npart(nn)     ! element and node partition vectors

    integer(idx_t) :: opts(0:METIS_NOPTIONS-1), ios, objval, nparts

    write(*,'(A)') "TEST METIS_PartMeshNodal2"

    !  1---2---5
    !  | 1 | 2 |
    !  4---3---6
    !  | 4 | 3 |
    !  9---8---7

    eptr = [1,5,9,13,17]
    eind = [1,2,3,4,&
            2,5,6,3,&
            3,6,7,8,&
            4,3,8,9]

    ios = METIS_SetDefaultOptions(opts)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
        error stop 1
    end if
    opts(METIS_OPTION_NUMBERING) = 1    ! Fortran-style numbering
    opts(METIS_OPTION_CONTIG) = 1       ! Force contigous partitions
    ! call print_metis_options(opts)
    
    nparts = 2
    ios = METIS_PartMeshNodal(ne,nn,eptr,eind,nparts=nparts,options=opts, &
            objval=objval,epart=epart,npart=npart)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_PartMeshNodal failed with error: ", ios
        error stop 1
    end if

    write(*,'(A,I0)') "objval = ", objval
    write(*,'(A,*(I1,:,1X))') "epart = ", epart
    write(*,'(A,*(I1,:,1X))') "npart = ", npart

end program