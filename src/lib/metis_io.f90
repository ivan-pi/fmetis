! metis_io.f90 -- Fortran METIS Interface
!
! Copyright (C) 2018 Ivan Pribec <ivan.pribec@gmail.com>
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

module metis_io

    use metis_interface
    use metis_enum

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

            do i = 1, size(xadj)-1
                write(unit,*) (adjncy(j), j = xadj(i)+1, xadj(i+1))
            end do
        else
            write(unit,*) size(xadj)-1,size(adjncy)/2
            do i = 1, size(xadj)-1
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
        integer :: unit, ncol, ios, i, rowcol, j
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

    subroutine print_metis_options(opts,unit)
        use iso_fortran_env, only: output_unit
        integer, intent(in) :: opts(0:)
        integer, intent(in), optional :: unit
        integer :: i, unit_

        unit_ = output_unit ! standard output
        if (present(unit)) unit_ = unit

        do i = 0, METIS_NOPTIONS-1
            write(unit_,'("Option ",I2,":",I3)') i, opts(i)
        end do
    end subroutine


    subroutine ForMETIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy,stat)
        use iso_c_binding, only : c_int, c_ptr, c_f_pointer
        integer, intent(in) :: ne
        integer, intent(in) :: nn
        integer, intent(in) :: eptr(ne+1)
        integer, intent(in) :: eind(:)
        integer, intent(in) :: numflag
        integer, intent(out), allocatable :: xadj(:)
        integer, intent(out), allocatable :: adjncy(:)
        integer, intent(out) :: stat

        type(c_ptr) :: xadj_, adjncy_
        
        integer(c_int), pointer :: fxadj(:) => null(), fadjncy(:) => null()

        stat = METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj_,adjncy_)
        if (stat < 0) return

        call c_f_pointer(xadj_,fxadj,shape=[nn+1])
        
        select case(numflag)
        case(0)
            call c_f_pointer(adjncy_,fadjncy,shape=[fxadj(nn+1)])
        case(1)
            call c_f_pointer(adjncy_,fadjncy,shape=[fxadj(nn+1)-1])
        end select

        ! xadj => fxadj
        allocate(xadj(size(fxadj)))
        xadj = fxadj
        allocate(adjncy(size(fadjncy)),stat=stat)
        adjncy = fadjncy

        stat = METIS_Free(xadj_); 
        ! if (stat < 0) return
        stat = METIS_Free(adjncy_)
        ! if (stat < 0) return
    end subroutine

end module