!>  Partition a graph specified by a matrix
!
!   Source: http://people.eecs.berkeley.edu/~demmel/cs267/lecture18/lecture18.html
!
program test_PartGraphKway

    use metis_interface, only: idx_t, METIS_PartGraphKway, METIS_SetDefaultOptions, &
        METIS_NOPTIONS, METIS_OPTION_NUMBERING, METIS_OK
    implicit none

    integer(idx_t), parameter :: npart = 3  ! number of partitions
    integer(idx_t), parameter :: n = 8      ! number of nodes
    integer(idx_t), parameter :: m = 10     ! number of edges

    integer(idx_t) :: xadj(n+1), adjncy(2*m) ! graph adjacency structure
    integer(idx_t) :: part(n)
    integer(idx_t) :: options(0:METIS_NOPTIONS-1), ios, objval

    write(*,'(A)') "TEST METIS_PartGraphKway"

    !    1  2  3  4  5  6  7  8
    !  | a  a  a     a          | 1
    !  | a  a     a  a          | 2
    !  | b     b        b       | 3
    !  |    b     b     b       | 4
    !  | a  a           a     a | 5
    !  |       b  b  b          | 6
    !  |                   c  c | 7
    !  |             c     c  c | 8

    xadj = [1,4,7,9,11,15,18,19,21]
    adjncy = [2,3,5,&
              1,4,5,&
              1,6,&
              2,6,&
              1,2,6,8,&
              3,4,5,&
              8,&
              5,7]

    ios = METIS_SetDefaultOptions(options)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
        error stop 1
    end if
    options(METIS_OPTION_NUMBERING) = 1 ! Fortran-style numbering

    ios = METIS_PartGraphKway(n,ncon=1,xadj=xadj,adjncy=adjncy,nparts=npart,options=options,objval=objval,part=part)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_PartGraphKway failed with error: ", ios
        error stop 1
    end if

    write(*,'(A,I0)') "objval = ", objval
    write(*,'(A,*(I1,:,1X))') "part = ", part

end program