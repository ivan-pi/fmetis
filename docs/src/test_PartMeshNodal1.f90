! https://stackoverflow.com/questions/20006253/using-metis-libraries-in-fortran-code-the-basics
! http://glaros.dtc.umn.edu/gkhome/node/852

program test_PartMeshNodal1

    use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
    use metis_interface, only: METIS_SetDefaultOptions, METIS_PartMeshNodal, &
        METIS_MeshToNodal, METIS_Free, METIS_NOPTIONS
    use metis_enum, only: METIS_OK, METIS_OPTION_NUMBERING
    implicit none

    integer(c_int), parameter :: ne = 3 ! number of elements
    integer(c_int), parameter :: nn = 8 ! number of nodes
    integer(c_int), parameter :: npel = 4 ! nodes per element

    integer(c_int) :: eptr(ne+1)
    integer(c_int) :: eind(ne*npel)
    integer(c_int) :: epart(ne), npart(nn)
    integer(c_int) :: options(0:METIS_NOPTIONS-1)
    integer(c_int) :: ios, objval

    type(c_ptr) :: c_xadj, c_adjncy
    integer(c_int), pointer :: xadj(:) => null(), adjncy(:) => null()

    write(*,'(A)') "TEST METIS_PartMeshNodal 1"

    ! 0---1---4---6
    ! | 0 | 1 | 2 |
    ! 3---2---5---7

    eptr = [0,4,8,12]
    eind = [0,1,2,3,1,4,5,2,4,6,7,5]

    ios = METIS_SetDefaultOptions(options)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
        error stop 1
    end if
    options(METIS_OPTION_NUMBERING) = 0 ! C-style numbering

    ios = METIS_PartMeshNodal(ne,nn,eptr,eind,nparts=2,options=options,&
        objval=objval,epart=epart,npart=npart)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_PartMeshNodal failed with error: ", ios
        error stop 1
    end if

    write(*,'(A,I0)') "objval = ", objval
    write(*,'(A,*(I1,:,1X))') "epart = ", epart
    write(*,'(A,*(I1,:,1X))') "npart = ", npart

    ios = METIS_MeshToNodal(ne,nn,eptr,eind,numflag=0,xadj=c_xadj,adjncy=c_adjncy)
    if (ios /= METIS_OK) then
        write(*,*) "METIS_MeshToNodal failed with error: ", ios
        error stop 1
    end if

    call c_f_pointer(c_xadj,xadj,shape=[nn+1]) ! size of adjacency list is one more than number of nodes
    call c_f_pointer(c_adjncy,adjncy,shape=[xadj(nn+1)]) ! size of edge list is in the last element of xadj

    write(*,'(A,*(I0,:,1X))') "xadj = ", xadj
    write(*,'(A,*(I1,:,1X))') "adjncy = ", adjncy

    ! call write_graph("test1.graph",xadj,adjncy,1)

    ios = METIS_Free(c_xadj); xadj => null()
    if (ios /= METIS_OK) then
        write(*,*) "METIS_Free failed with error: ", ios
        error stop 1
    end if
    ios = METIS_Free(c_adjncy); adjncy => null()
    if (ios /= METIS_OK) then
        write(*,*) "METIS_Free failed with error: ", ios
        error stop 1
    end if

end program