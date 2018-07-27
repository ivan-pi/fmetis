module metis_tests

    use metis_interface
    use metis_io

    implicit none
    private

    public :: test1, test2, test3, test4

contains
    ! https://stackoverflow.com/questions/20006253/using-metis-libraries-in-fortran-code-the-basics
    ! http://glaros.dtc.umn.edu/gkhome/node/852
    subroutine test1()
        use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer

        integer, parameter :: nel = 3 ! number of elements
        integer, parameter :: nnds = 8 ! number of nodes
        integer, parameter :: npel = 4 ! nodes per element

        integer :: eptr(nel+1)
        integer :: eind(nel*npel)
        integer :: epart(nel), npart(nnds)
        integer(c_int) :: options(0:METIS_NOPTIONS-1)
        integer :: ios, objval

        type(c_ptr) :: xadj, adjncy
        
        integer(c_int), pointer :: fxadj(:) => null(), fadjncy(:) => null()
        integer, allocatable :: new_xadj(:), new_adjncy(:)

        print *, "TEST 1"

        ! 0---1---4---6
        ! : 0 : 1 : 2 :
        ! 3---2---5---7

        eptr = [0,4,8,12]
        eind = [0,1,2,3,1,4,5,2,4,6,7,5]  ! Element 1 has nodes 0 1 2 3
                                 ! Element 2 has nodes 1 4 5 2

        ios = METIS_SetDefaultOptions(options)
        options(17) = 0

        ios = METIS_PartMeshNodal(nel,nnds,eptr,eind,nparts=2,options=options,&
            objval=objval,epart=epart,npart=npart)

        print *, "ios = ", ios
        print *, "objval = ", objval
        print *, "npart = ", npart
        print *, "epart = ", epart


        call ForMETIS_MeshToNodal(nel,nnds,eptr,eind,0,new_xadj,new_adjncy,ios)
        print *, "new_xadj = ", new_xadj
        print *, "new_adjncy = ", new_adjncy

        ios = METIS_MeshToNodal(nel,nnds,eptr,eind,0,xadj,adjncy)
        ! print *, "xadj = ", xadj
        ! print *, "adjncy = ", adjncy

        call c_f_pointer(xadj,fxadj,shape=[nnds+1])
        call c_f_pointer(adjncy,fadjncy,shape=[fxadj(nnds+1)])

        print *, "size=",size(fxadj),"fxadj = ", fxadj
        print *, "fadjncy = ", fadjncy

        call write_graph("test1.graph",fxadj,fadjncy,1)

        ios = METIS_Free(xadj)
        ! print *, "xadj = ", xadj

        ! print *, "fxadj = ", fxadj

        ! ios = METIS_Free(adjncy)
    end subroutine

    ! https://stackoverflow.com/questions/8155160/metis-with-fortran
    ! http://glaros.dtc.umn.edu/gkhome/node/799
    subroutine test2()
        use, intrinsic :: iso_c_binding, only : c_int

        integer, parameter :: nvtxs = 15
        integer, parameter :: nedgs = 22

        integer :: xadj(nvtxs+1),adjncy(2*nedgs)
        integer :: objval, part(nvtxs), ios

        print *, "TEST 2"

        xadj = [0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44]
        adjncy = [1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13]

        ! options(18) = 0 ! C-style indexing

        ios = METIS_PartGraphRecursive(nvtxs,ncon=1,xadj=xadj,adjncy=adjncy,nparts=2,objval=objval,part=part)

        print *, "ios = ", ios
        print *, "objval = ", objval
        print *, "part = ", part
    end subroutine

    ! http://comp.lang.fortran.narkive.com/uFmDM7Bo/how-to-call-a-metis-subroutine-from-my-fortran-code
    subroutine test3()

        !  1---2---3---4---5
        !  |   |   |   |   |
        !  6---7---8---9---10
        !  |   |   |   |   |
        !  11--12--13--14--15

        integer, parameter :: n = 15
        integer, parameter :: m = 22

        integer :: xadj(n+1), adjncy(2*m)
        integer :: perm(n), iperm(n)

        integer :: options(0:40), ios

        print *, "TEST 3"

        xadj = [1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45]
        adjncy = [2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6, &
                  8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,7,11,13, &
                  8,12,14,9,13,15,10,14]

        ios = METIS_SetDefaultOptions(options)
        options(17) = 1

        ios = METIS_NodeND(n,xadj,adjncy,options=options,perm=perm,iperm=iperm)

        print *, "ios = ", ios
        print *, "perm = ", perm
        print *, "iperm = ", iperm

    end subroutine

    ! https://www.cfd-online.com/Forums/main/112366-using-metis-functions-fortran.html
    subroutine test4

        !  1---2---5
        !  | 1 | 2 |
        !  4---3---6
        !  | 4 | 3 |
        !  9---8---7

        integer, parameter :: ne = 4
        integer, parameter :: nn = 9

        integer :: eptr(ne+1), eind(4*ne)
        integer :: epart(ne), npart(nn)

        integer :: opts(0:39), ios, objval

        print *, "TEST 4"

        eptr = [1,5,9,13,17]
        eind = [1,2,3,4,2,5,6,3,3,6,7,8,4,3,8,9]

        ios = METIS_SetDefaultOptions(opts)
        opts(17) = 1
        opts(METIS_OPTION_CONTIG) = 1
        call print_options(opts)
        
        ios = METIS_PartMeshDual(ne,nn,eptr,eind,ncommon=2,nparts=2,options=opts, &
            objval=objval,epart=epart,npart=npart)

        ! ios = METIS_PartMeshNodal(ne,nn,eptr,eind,nparts=2,options=opts, &
            ! objval=objval,epart=epart,npart=npart)

        print *, "ios = ", ios
        print *, "objval = ", objval
        print *, "epart = ", epart
        print *, "npart = ", npart

    end subroutine

end module

program test_fmetis

    use iso_c_binding
    use metis_interface
    use metis_enum
    use metis_io
    use metis_tests, only : test1, test2, test3, test4

    type(MetisGraph) :: mgraph
    integer, allocatable :: xadj(:), adjncy(:), part(:)

    integer(c_int) :: options(40), ios, ncon, objval

    xadj = [0,3,6,10,14,17,20,22]

    adjncy = [5,3,2, 1,3,4, 5,4,2,1 , 2,3,6,7, 1,3,6, 5,4,7, 6,4] - 1

    print *, "n = ", size(xadj) -1
    print *, "2*m = ", size(adjncy)

    call write_graph("example.graph",xadj,adjncy,0)

    ! call load_graph(mgraph,"../graphs/4elt.graph",numflag=0)
    call load_graph(mgraph,"example.graph",numflag=0)

    print *, mgraph%nvtxs
    print *, mgraph%nedgs
    print *, mgraph%xadj
    print *, mgraph%adjncy
    ! call write_graph("print.graph",)

    ios = METIS_SetDefaultOptions(options)
    options(18) = 0
    
    print *, "ios = ", ios

    allocate(part(mgraph%nvtxs))

    ncon = 1
    ios = METIS_PartGraphKway(mgraph%nvtxs,ncon,mgraph%xadj,mgraph%adjncy,&
        nparts=2,objval=objval,part=part,options=options)

    print *, "ios = ", ios
    print *, "objval = ", objval
    print *, "part = ", part


    write(*,*)
    call test1()

    print *, " " 
    call test2()

    print *, " "
    call test3()

    print *, " "
    call test4()
end program