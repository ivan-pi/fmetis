module fmetis_interface

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
    public :: METIS_NOPTIONS

    integer, parameter :: METIS_NOPTIONS = 40

    ! http://glaros.dtc.umn.edu/gkhome/node/877

    interface

        integer(c_int) function METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,&
                                                         vwgt,vsize,adjwgt,nparts,tpwgts,&
                                                         ubvec,options,objval,part) bind(C,name="METIS_PartGraphRecursive")
            use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr
            import METIS_NOPTIONS
            integer(c_int), intent(in) :: nvtxs
                !! The number of vertices in the graph.
            integer(c_int), intent(in) :: ncon
                !! The number of balancing constraints. It should be atleast 1.
            integer(c_int), intent(in), dimension(*) :: xadj, adjncy
                !! The adjacency structure of the graph.
            integer(c_int), intent(in), dimension(*), optional :: vwgt ! NULL
                !! The weights of the vertices.
            integer(c_int), intent(in), dimension(*), optional :: vsize ! NULL
                !! The size of the vertices for computing the total communication volume.
            integer(c_int), intent(in), dimension(*), optional :: adjwgt ! NULL
                !! The weights of the edges.
            integer(c_int) , intent(in) :: nparts
                !! The number of parts to partition the graph.
            real(c_double), dimension(*), optional :: tpwgts
                !! The target partition weights.
            real(c_double), dimension(*), optional :: ubvec
                !! Array for specifying the allowed load imbalance
            integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            integer(c_int), intent(out) :: objval
            integer(c_int), intent(out) :: part(nvtxs+1)
        end function

        integer(c_int) function METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,&
                                                         vwgt,vsize,adjwgt,nparts,tpwgts,&
                                                         ubvec,options,objval,part) bind(C,name="METIS_PartGraphKway")
            use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr
            import METIS_NOPTIONS
            integer(c_int), intent(in) :: nvtxs
                !! The number of vertices in the graph.
            integer(c_int), intent(in) :: ncon
                !! The number of balancing constraints. It should be atleast 1.
            integer(c_int), intent(in), dimension(*) :: xadj, adjncy
                !! The adjacency structure of the graph.
            integer(c_int), intent(in), dimension(*), optional :: vwgt ! NULL
                !! The weights of the vertices.
            integer(c_int), intent(in), dimension(*), optional :: vsize ! NULL
                !! The size of the vertices for computing the total communication volume.
            integer(c_int), intent(in), dimension(*), optional :: adjwgt ! NULL
                !! The weights of the edges.
            integer(c_int), intent(in) :: nparts
                !! The number of parts to partition the graph.
            real(c_double), dimension(*), optional :: tpwgts
                !! The target partition weights.
            real(c_double), dimension(*), optional :: ubvec
                !! Array for specifying the allowed load imbalance
            integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            integer(c_int), intent(out) :: objval
            integer(c_int), intent(out) :: part(nvtxs)
        end function

        integer(c_int) function METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon, &
            nparts,tpwgts,options,objval,epart,npart) bind(C,name="METIS_PartMeshDual")
            use, intrinsic :: iso_c_binding, only : c_int
            import METIS_NOPTIONS
            integer(c_int), intent(in) :: ne
                !! The number of elements in the mesh.
            integer(c_int), intent(in) :: nn
                !! The number of nodes in the mesh.
            integer(c_int), intent(in), dimension(*) :: eptr,eind
                !! The pair of arrays storing the mesh.
            integer(c_int), intent(in), dimension(*), optional :: vwgt
            integer(c_int), intent(in), dimension(*), optional :: vsize
            integer(c_int), intent(in) :: ncommon
            integer(c_int), intent(in) :: nparts
            integer(c_int), dimension(*), optional :: tpwgts
            integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            integer(c_int), intent(out) :: objval
            integer(c_int), intent(out) :: epart(ne)
            integer(c_int), intent(out) :: npart(nn)
        end function


        integer(c_int) function METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsize, &
            nparts,tpwgts,options,objval,epart,npart) bind(C,name="METIS_PartMeshNodal")
            use, intrinsic :: iso_c_binding, only : c_int
            import METIS_NOPTIONS
            integer(c_int), intent(in) :: ne
                !! The number of elements in the mesh.
            integer(c_int), intent(in) :: nn
                !! The number of nodes in the mesh.
            integer(c_int), intent(in), dimension(*) :: eptr,eind
                !! The pair of arrays storing the mesh.
            integer(c_int), intent(in), dimension(*), optional :: vwgt
            integer(c_int), intent(in), dimension(*), optional :: vsize
            integer(c_int), intent(in) :: nparts
            integer(c_int), dimension(*), optional :: tpwgts
            integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            integer(c_int), intent(out) :: objval
            integer(c_int), intent(out) :: epart(ne)
            integer(c_int), intent(out) :: npart(nn)
        end function

        integer(c_int) function METIS_NodeND(nvtxs,xadj,adjncy,vwgt,options,perm,iperm) bind(C,name="METIS_NodeND")
            use iso_c_binding, only : c_int, c_ptr
            import METIS_NOPTIONS
            integer(c_int), intent(in) :: nvtxs
            integer(c_int), intent(in), dimension(*) :: xadj, adjncy
            integer(c_int), intent(in), optional :: vwgt(nvtxs)
            integer(c_int), intent(in), optional :: options(METIS_NOPTIONS)
            integer(c_int), intent(out) :: perm(nvtxs), iperm(nvtxs)
        end function


        integer(c_int) function METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy) bind(C,name="METIS_MeshToDual")
            use, intrinsic :: iso_c_binding, only : c_int, c_ptr
            integer(c_int), intent(in) :: ne, nn
            integer(c_int), intent(in), dimension(*) :: eptr, eind
            integer(c_int), intent(in) :: ncommon
            integer(c_int), intent(in) :: numflag
                !! Used to indicate which numbering scheme is used for `eptr` and `eind`.
                !! The possible values are: 
                !!  * 0 - C-style numbering is assumed that starts from 0
                !!  * 1 - Fortran-style numbering is assumed that starts from 1
            type(c_ptr), intent(out) :: xadj, adjncy
            ! integer(c_int), allocatable :: adjncy(:)
        end function

        integer(c_int) function METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy) bind(C,name="METIS_MeshToNodal")
            use, intrinsic :: iso_c_binding, only : c_int, c_ptr
            integer(c_int), intent(in) :: ne, nn
            integer(c_int), intent(in), dimension(*) :: eptr, eind
            integer(c_int), intent(in) :: numflag
                !! Used to indicate which numbering scheme is used for `eptr` and `eind`.
                !! The possible values are: 
                !!  * 0 - C-style numbering is assumed that starts from 0
                !!  * 1 - Fortran-style numbering is assumed that starts from 1
            type(c_ptr), intent(out) :: xadj, adjncy
            ! integer(c_int), allocatable :: adjncy(:)
        end function

        integer(c_int) function METIS_SetDefaultOptions(options) bind(C,name="METIS_SetDefaultOptions")
            use, intrinsic :: iso_c_binding, only : c_int, c_ptr
            import METIS_NOPTIONS
            integer(c_int), intent(out) :: options(METIS_NOPTIONS)
        end function

        integer(c_int) function METIS_Free(ptr) bind(c)
            use iso_c_binding, only : c_int, c_ptr
            type(c_ptr), value :: ptr
        end function

    end interface ! METIS


end module

module fmetis_wrapper

    use fmetis_interface
    implicit none

    public

    type :: MetisGraph
        integer :: nvtxs, nedgs
        integer, allocatable :: xadj(:), adjncy(:)
    end type

contains

    subroutine write_graph(fname,xadj,adjncy,numflag)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: xadj(:)
        integer, intent(in) :: adjncy(:)
        integer, intent(in) :: numflag

        integer :: unit, i, j

        print *, "[write_graph] adjncy = ", adjncy

        open(newunit=unit,file=fname)

        if (numflag == 0) then
            write(unit,*) size(xadj)-1, size(adjncy)/2

            do i = 1, size(xadj)
                write(unit,*) (adjncy(j), j = xadj(i)+1, xadj(i+1))
            end do
        else
            write(unit,*) size(xadj)-1,size(adjncy)/2
            do i = 1, size(xadj)
                write(unit,*) (adjncy(j), j = xadj(i), xadj(i+1)-1)
            end do
        end if

        close(unit)
    end subroutine

    logical function whitechar(char) ! white character
    ! returns .true. if char is space (32) or tab (9), .false. otherwise
    character, intent(in) :: char
    if (iachar(char) == 32 .or. iachar(char) == 9) then
        whitechar = .true.
    else
        whitechar = .false.
    end if
    end function

    integer function count_columns(unit,stat) result(ncol)
        integer, intent(in) :: unit
        integer, intent(out) :: stat
        
        character(len=1) :: c
        logical :: lastwhite

        ncol = 0
        lastwhite = .true.
        do
            read(unit,'(a)',advance='no',iostat=stat) c
            if (stat /= 0) exit
            if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
            lastwhite = whitechar(c)
        end do
    end function

    subroutine load_graph(this,fname,numflag)
        type(MetisGraph), intent(out) :: this
        character(len=*), intent(in) :: fname
        integer, intent(in) :: numflag

        character(len=1) :: c
        integer :: unit, ncol, ios, i, rowcol, j, k
        logical :: lastwhite

        character(len=3) :: fmt
        integer :: n, m, ncon, ifmt

        integer, allocatable :: xadj(:), adjncy(:)
        integer, allocatable :: vsize(:), vwgt(:), adjwgt(:)
        integer, allocatable :: tmp(:)

        open(newunit=unit, file=fname, status='old')

        ! determine number of columns in first row of file
        ncol = 0
        lastwhite = .true.
        do
            read(unit,'(a)',advance='no',iostat=ios) c
            if (ios /= 0) exit
            if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
            lastwhite = whitechar(c)
        end do

        ! back to start of file
        rewind(unit)

        ! read header line
        select case(ncol)
        case(2)
            fmt = "000"
            ncon = 0
            read(unit,*,iostat=ios) n, m
        case(3)
            ncon = 0
            read(unit,*,iostat=ios) n, m, fmt
        case(4)
            read(unit,*,iostat=ios) n, m, fmt, ncon
        case default
            write(*,*) "[load_graph]: incorrect file"
            stop
        end select
        read(fmt,'(I3)',iostat=ios) ifmt

        print *, n, m, fmt, ncon
        print *, "ifmt = ", ifmt

        ! allocate vertex adjacency structures
        allocate(xadj(n+1))
        allocate(adjncy(2*m))

        ! perform bit-tests to check for weights and sizes
        if (btest(ifmt,0)) allocate(adjwgt(2*m))
        if (btest(ifmt,1)) then
            if (ncon == 0) then
                ncol = 4 ! implicit fourth column with 1 constraint
                ncon = 1 ! case for fmt = "*1*" and ncon is not specified
            end if
            if (ncon > 0) then
                allocate(vwgt(n*ncon))
            else
                print *, "[load_graph] ncon has not been specified"
                stop
            end if
        end if
        if (btest(ifmt,2)) allocate(vsize(n))
        
        print *, "allocated vertex sizes", allocated(vsize)
        print *, "allocated vertex weights", allocated(vwgt)
        print *, "allocated edge weights", allocated(adjwgt)

        xadj(1) = 0
        select case(ncol)
        case (2) ! only connectivity, fmt = "000"
            do i = 1, n
                rowcol = count_columns(unit,stat=ios)
                ! print *, "rowcol = ", rowcol
                xadj(i+1) = xadj(i) + rowcol
                backspace(unit,iostat=ios)
                read(unit,*) adjncy(xadj(i)+1:xadj(i+1))
                ! print *, (adjncy(j),j=xadj(i)+1, xadj(i+1))
            end do
        case (3) ! edge weights or vertex sizes
            select case(fmt)
            case("001")
                do i = 1, n
                    rowcol = count_columns(unit,stat=ios)
                    backspace(unit,iostat=ios)
                    ! print *, "rowcol = ", rowcol
                    xadj(i+1) = xadj(i) + rowcol/2
                    
                    if (allocated(tmp)) deallocate(tmp)
                    allocate(tmp(rowcol))
                    
                    read(unit,*) tmp
                    
                    adjncy(xadj(i)+1:xadj(i+1)) = tmp(1::2)
                    adjwgt(xadj(i)+1:xadj(i+1)) = tmp(2::2)
                end do
            case("101")
                do i = 1, n
                    rowcol = count_columns(unit,stat=ios) - 1
                    backspace(unit,iostat=ios)
                    ! print *, "rowcol = ", rowcol
                    xadj(i+1) = xadj(i) + rowcol/2
                    
                    if (allocated(tmp)) deallocate(tmp)
                    allocate(tmp(rowcol))
                    
                    read(unit,*) vsize(i), tmp
                    
                    adjncy(xadj(i)+1:xadj(i+1)) = tmp(1::2)
                    adjwgt(xadj(i)+1:xadj(i+1)) = tmp(2::2)
                end do
            case("100")
                do i = 1, n
                    rowcol = count_columns(unit,stat=ios) - 1
                    ! print *, "rowcol = ", rowcol
                    xadj(i+1) = xadj(i) + rowcol
                    backspace(unit,iostat=ios)
                    read(unit,*) vsize(i), (adjncy(j),j=xadj(i)+1, xadj(i+1))
                    ! print *, (adjncy(j),j=xadj(i)+1, xadj(i+1))
                end do
            case default
                print *, "should not be here"
                stop
            end select
        case (4) ! vertex weights and constraint
            select case(fmt)
            case('010')
                do i = 1, n
                    rowcol = count_columns(unit,stat=ios) - ncon
                    print *, ncon, rowcol
                    
                    xadj(i+1) = xadj(i) + rowcol
                    backspace(unit,iostat=ios)
                    read(unit,*) vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), adjncy(xadj(i)+1:xadj(i+1))
                end do
            case('011')
                do i = 1, n
                    rowcol = (count_columns(unit,stat=ios) - ncon)/2
                    print *, ncon, rowcol
                    
                    backspace(unit,iostat=ios)
                    xadj(i+1) = xadj(i) + rowcol
                    
                    if (allocated(tmp)) deallocate(tmp)
                    allocate(tmp(2*rowcol))
                    
                    read(unit,*) vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), tmp
                    adjncy(xadj(i)+1:xadj(i+1)) = tmp(1::2)
                    adjwgt(xadj(i)+1:xadj(i+1)) = tmp(2::2)
                end do
            case('110')
                do i = 1, n
                    rowcol = count_columns(unit,stat=ios) - 1 - ncon
                    print *, ncon, rowcol
                    
                    xadj(i+1) = xadj(i) + rowcol
                    backspace(unit,iostat=ios)
                    read(unit,*) vsize(i), vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), adjncy(xadj(i)+1:xadj(i+1))
                end do

            case('111')
                do i = 1, n
                    rowcol = (count_columns(unit,stat=ios) - 1 - ncon)/2
                    print *, ncon, rowcol
                    
                    backspace(unit,iostat=ios)
                    xadj(i+1) = xadj(i) + rowcol
                    
                    if (allocated(tmp)) deallocate(tmp)
                    allocate(tmp(2*rowcol))
                    
                    read(unit,*) vsize(i), vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), tmp
                    adjncy(xadj(i)+1:xadj(i+1)) = tmp(1::2)
                    adjwgt(xadj(i)+1:xadj(i+1)) = tmp(2::2)
                end do
            case default
                print *, "should not be here"
                stop
            end select ! fmt
        case default 
            print *, "[load_graph] error"
            stop
        end select ! ncol

        this%nvtxs = n
        this%nedgs = m
        this%xadj = xadj    

        if (numflag == 0) then
            this%adjncy = adjncy
        else
            this%adjncy = adjncy + 1
        end if

        ! print *, "xadj = ", xadj
        ! print *, "adjncy = ",adjncy
        ! print *, "edge weights = ", adjwgt
        ! print *, "vsize = ", vsize
        ! print *, "vertex weights = ", vwgt
        close(unit)

    end subroutine

    subroutine print_options(opts)
        integer, intent(in) :: opts(:)
        integer :: i

        do i = lbound(opts,1), ubound(opts,1)
            print *, "option ",i, opts(i)
        end do
    end subroutine

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
        integer(c_int) :: options(METIS_NOPTIONS)
        integer :: ios, objval

        type(c_ptr) :: xadj, adjncy
        
        integer(c_int), pointer :: fxadj(:) => null(), fadjncy(:) => null()

        print *, "TEST 1"

        ! 0---1---4---6
        ! : 0 : 1 : 2 :
        ! 3---2---5---7

        eptr = [0,4,8,12] + 1
        eind = [0,1,2,3,1,4,5,2,4,6,7,5] + 1  ! Element 1 has nodes 0 1 2 3
                                 ! Element 2 has nodes 1 4 5 2

        ios = METIS_SetDefaultOptions(options)
        options(18) = 1

        ios = METIS_PartMeshNodal(nel,nnds,eptr,eind,nparts=2,options=options,&
            objval=objval,epart=epart,npart=npart)

        print *, "ios = ", ios
        print *, "objval = ", objval
        print *, "npart = ", npart
        print *, "epart = ", epart

        ios = METIS_MeshToNodal(nel,nnds,eptr,eind,1,xadj,adjncy)
        print *, "xadj = ", xadj
        print *, "adjncy = ", adjncy

        call c_f_pointer(xadj,fxadj,shape=[nnds+1])
        call c_f_pointer(adjncy,fadjncy,shape=[fxadj(nnds+1)])

        print *, "fxadj = ", fxadj
        print *, "fadjncy = ", fadjncy

        call write_graph("test1.graph",fxadj,fadjncy,1)

        ios = METIS_Free(xadj)
        print *, "xadj = ", xadj

        print *, "fxadj = ", fxadj

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

        integer :: options(40), ios

        print *, "TEST 3"

        xadj = [1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45]
        adjncy = [2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6, &
                  8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,7,11,13, &
                  8,12,14,9,13,15,10,14]

        ios = METIS_SetDefaultOptions(options)
        options(18) = 1

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

        integer :: opts(40), ios, objval

        print *, "TEST 4"

        eptr = [1,5,9,13,17]
        eind = [1,2,3,4,2,5,6,3,3,6,7,8,4,3,8,9]

        ios = METIS_SetDefaultOptions(opts)
        opts(18) = 1

        ! ios = METIS_PartMeshDual(ne,nn,eptr,eind,ncommon=2,nparts=2,options=opts, &
            ! objval=objval,epart=epart,npart=npart)

        ios = METIS_PartMeshNodal(ne,nn,eptr,eind,nparts=2,options=opts, &
            objval=objval,epart=epart,npart=npart)

        print *, "ios = ", ios
        print *, "objval = ", objval
        print *, "epart = ", epart
        print *, "npart = ", npart

    end subroutine

end module

program test_fmetis

    use iso_c_binding
    use fmetis_interface
    use fmetis_wrapper

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


    print *,
    call test1

    print *, 
    call test2

    print *,
    call test3

    print *,
    call test4
end program

