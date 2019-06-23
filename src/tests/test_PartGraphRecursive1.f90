! https://stackoverflow.com/questions/8155160/metis-with-fortran
! http://glaros.dtc.umn.edu/gkhome/node/799

program test_PartGraphRecursive1

    use metis_interface, only: idx_t, METIS_PartGraphRecursive, METIS_OK
    implicit none
    
    integer(idx_t), parameter :: nvtxs = 15     ! number of vertices
    integer(idx_t), parameter :: nedgs = 22     ! number of edges

    integer(idx_t) :: xadj(nvtxs+1),adjncy(2*nedgs) ! adjacency arrays
    integer(idx_t) :: part(nvtxs)                   ! partition vector
    integer(idx_t) :: objval, ios, ncon, nparts

    write(*,'(A)') "TEST METIS_PartGraphRecursive 1"

    !  0---1---2---3---4
    !  |   |   |   |   |
    !  5---6---7---8---9
    !  |   |   |   |   |
    ! 10--11--12--13--14

    ! Note that we are using C-style numbering!

    xadj = [0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44]
    adjncy = [1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13]

    ncon = 1
    nparts = 2
    ios = METIS_PartGraphRecursive(nvtxs,ncon,xadj=xadj,adjncy=adjncy,nparts=nparts,objval=objval,part=part)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_PartGraphRecursive failed with error: ", ios
        error stop 1
    end if

    write(*,'(A,I0)') "objval = ", objval
    write(*,'(A,*(I1,:,1X))') "part = ", part

end program