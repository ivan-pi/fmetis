!>  Example of computing a fill-reducing ordering of a sparse matrix
!
!   Source: http://comp.lang.fortran.narkive.com/uFmDM7Bo/how-to-call-a-metis-subroutine-from-my-fortran-code
!
program test_NodeND

    use iso_c_binding, only: c_int
    use metis_interface, only: METIS_SetDefaultOptions, METIS_NodeND, METIS_NOPTIONS
    use metis_enum, only: METIS_OK, METIS_OPTION_NUMBERING
    implicit none

    integer(c_int), parameter :: n = 15 ! number of vertices
    integer(c_int), parameter :: m = 22 ! number of edges

    integer(c_int) :: xadj(n+1), adjncy(2*m) ! graph adjacency structure
    integer(c_int) :: perm(n), iperm(n) ! fill-reducing permutation and inverse permutation

    integer(c_int) :: options(0:METIS_NOPTIONS-1), ios

    write(*,'(A)') "TEST METIS_NodeND"

    !  1---2---3---4---5
    !  |   |   |   |   |
    !  6---7---8---9---10
    !  |   |   |   |   |
    !  11--12--13--14--15

    xadj = [1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45]
    adjncy = [2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6, &
              8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,7,11,13, &
              8,12,14,9,13,15,10,14]

    ios = METIS_SetDefaultOptions(options)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
        error stop 1
    end if
    options(METIS_OPTION_NUMBERING) = 1 ! Fortran-style numbering

    ios = METIS_NodeND(n,xadj,adjncy,options=options,perm=perm,iperm=iperm)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_NodeND failed with error: ", ios
        error stop 1
    end if

    write(*,'(A,*(I2,:,1X))') "perm  = ", perm
    write(*,'(A,*(I2,:,1X))') "iperm = ", iperm

end program