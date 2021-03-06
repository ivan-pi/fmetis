!>  Partitioning a graph with non-equal number of edges per node.
!
program test_PartGraphRecursive2

    use metis_interface, only: idx_t, METIS_PartGraphRecursive, METIS_SetDefaultOptions, &
        METIS_NOPTIONS, METIS_OPTION_NUMBERING, METIS_OK

    integer(idx_t) :: n, m
    integer(idx_t), allocatable :: xadj(:), adjncy(:), part(:)
    integer(idx_t) :: options(0:METIS_NOPTIONS-1), ios, ncon, objval, nparts

    write(*,'(A)') "TEST METIS_PartGraphRecursive 2"

    !   1---5
    !   |\ / \
    !   | 3   6      
    !   |/ \ / \
    !   2---4---7 

    xadj = [1,4,7,11,15,18,21,23]
    adjncy = [5,3,2,1,3,4,5,4,2,1,2,3,6,7,1,3,6,5,4,7,6,4]

    n = size(xadj) - 1
    m = size(adjncy)/2
    write(*,'(A,I0)') "n = ", n
    write(*,'(A,I0)') "m = ", m

    ios = METIS_SetDefaultOptions(options)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
        error stop 1
    end if
    options(METIS_OPTION_NUMBERING) = 1     ! Fortran-style numbering

    ncon = 1
    nparts = 2
    allocate(part(n))
    ios = METIS_PartGraphRecursive(n,ncon,xadj,adjncy,&
                nparts=nparts,objval=objval,part=part,options=options)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_PartGraphKway failed with error: ", ios
        error stop 1
    end if

    write(*,'(A,I0)') "objval = ", objval
    write(*,'(A,*(I1,:,1X))') "part = ", part

end program