! metis_oo_interface.f90 -- Fortran METIS Interface
!
! Copyright (C) 2018 Ivan Pribec <ivan.pribec@gmail.com>
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

module metis_oo_interface
    
    use metis_interface

    implicit none

    public :: graph_type
    public :: export_graph

    type :: graph_type
        integer :: nvxts !! Number of vertices.
        integer :: nedgs !! Number of edges.
        integer, pointer :: xadj(:) => null()
        integer, pointer :: adjncy(:) => null()
        integer :: numflag = 1 !! Numbering style.
        
        integer :: ncon
        integer, pointer :: vwgt(:) => null()
        integer, pointer :: adjwgt(:) => null()
        integer, pointer :: vsize(:) => null()
    end type

contains

    subroutine import_graph(fname,graph,numflag)
        character(len=*), intent(in) :: fname
        class(graph_type), intent(out) :: graph
        integer, intent(in), optional :: numflag

        integer :: unit,ios

        if (present(numflag)) graph%numflag = numflag
        
        open(newunit=unit,file=fname,status='old',iostat=ios)

        call read_graph(unit,graph%xadj,graph%adjncy,numflag=graph%numflag,vwgt=graph%vwgt, &
            adjwgt=graph%adjwgt,vsize=graph%vsize)

        print *, graph%adjncy

        graph%nvxts = size(graph%xadj)-1
        graph%nedgs = size(graph%adjncy)/2

        close(unit)
    end subroutine

    subroutine export_graph(fname,graph)
        character(len=*), intent(in) :: fname
        class(graph_type), intent(in) :: graph

        integer :: unit

        open(newunit=unit,file=fname)
        call write_graph(unit,graph%xadj,graph%adjncy,graph%numflag, &
            graph%vwgt,graph%adjwgt,graph%vsize)
        close(unit)
    end subroutine

    subroutine write_graph(unit,xadj,adjncy,numflag,vwgt,adjwgt,vsize)
        integer, intent(in) :: unit
        integer, intent(in) :: xadj(:)
        integer, intent(in) :: adjncy(:)
        integer, intent(in), optional :: numflag
        integer, intent(in), optional :: vwgt(:)
        integer, intent(in), optional :: adjwgt(:)
        integer, intent(in), optional :: vsize(:)

        integer :: nvxts, nedgs, ncon, numflag_, i, j, fmt
        character(len=3) :: cfmt
        character(len=11) :: fstring

        fstring = '(*(i0,:,x))' ! Format string for graph output

        numflag_ = 1 ! Assume Fortran numbering by default
        if (present(numflag)) numflag_ = numflag
        
        ! Get number of vertices and edges
        nvxts = size(xadj) - 1
        nedgs = size(adjncy)/2

        ! Format specifier
        fmt = 0
        if (present(adjwgt)) fmt = ibset(fmt,0)
        if (present(vwgt))   fmt = ibset(fmt,1)
        if (present(vsize))  fmt = ibset(fmt,2)

        ! Number of constraints
        ncon = 0
        if (btest(fmt,1)) ncon = size(vwgt)/nvxts

        ! Write header line
        if (fmt > 0) then
            ! Write fmt to character string
            write(cfmt,'(b3.3)') fmt 
            if (btest(fmt,1)) then
                if (ncon > 1) write(unit,'(i0,1x,i0,1x,a3,1x,i0)') nvxts, nedgs, cfmt, ncon
            else
                write(unit,'(i0,1x,i0,1x,a3)') nvxts, nedgs, cfmt
            end if
        else
            write(unit,'(i0,1x,i0)') nvxts, nedgs
        end if


        select case(fmt)
        case(b'000')
            do i = 1, nvxts
                ! v1 v2 v3 ...
                write(unit,fstring) (adjncy(j),j=xadj(i),xadj(i+1)-1)
            end do
        case(b'001') ! edge weights
            do i = 1, nvxts
                ! v1 e1 v2 e2 ...
                write(unit,fstring) (adjncy(j),adjwgt(j),j=xadj(i),xadj(i+1)-1)
            end do
        case(b'010') ! vertex weights
            do i = 1, nvxts
                ! w1 w2 ... wncon v1 v2 v3 ...
                write(unit,fstring) vwgt((i-1)*ncon+1:(i-1)*ncon+ncon),(adjncy(j),j=xadj(i),xadj(i+1)-1)
            end do
        case(b'100') ! vertex sizes
            do i = 1, nvxts
                write(unit,fstring) vsize(i), (adjncy(j),j=xadj(i),xadj(i+1)-1)
            end do
        case(b'011')
            do i = 1, nvxts
                write(unit,fstring) vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), (adjncy(j),adjwgt(j),j=xadj(i),xadj(i+1)-1)
            end do
        case(b'110')
            do i = 1, nvxts
                write(unit,fstring) vsize(i), vwgt((i-1)*ncon+1:(i-1)*ncon+ncon),(adjncy(j),j=xadj(i),xadj(i+1)-1)
            end do
        case(b'101')
            do i = 1, nvxts
                write(unit,fstring) vsize(i), (adjncy(j),adjwgt(j),j=xadj(i),xadj(i+1)-1)
            end do 
        case(b'111')
            do i = 1, nvxts
                write(unit,fstring) vsize(i), vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), (adjncy(j),adjwgt(j),j=xadj(i),xadj(i+1)-1)
            end do
        case default
            write(*,*) '[write_graph] Error occured'
        end select

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
        backspace(unit,iostat=stat)
    end function

    subroutine read_graph(unit,xadj,adjncy,numflag,vwgt,adjwgt,vsize)
        implicit none
        integer, intent(in) :: unit
        integer, intent(out), pointer :: xadj(:)
        integer, intent(out), pointer :: adjncy(:)
        integer, intent(in), optional :: numflag
        integer, pointer, optional :: vwgt(:)
        integer, intent(out), pointer, optional :: adjwgt(:)
        integer, intent(out), pointer, optional :: vsize(:)

        character(len=1) :: c
        integer :: ncol, ios, i, rowcol, j
        logical :: lastwhite

        character(len=3) :: cfmt
        integer :: nvtxs, nedgs, ncon, fmt, numflag_

        numflag_ = 1 ! Assume Fortran numbering by default
        if (present(numflag)) numflag_ = numflag

        ! Determine number of columns in header line
        ncol = 0
        lastwhite = .true.
        do
            read(unit,'(a)',advance='no',iostat=ios) c
            ! if (iachar(c) == 37) then
            !     read(unit,*) ! skipline
            !     print *, "Skipped line"
            !     cycle
            ! end if
            if (ios /= 0) exit
            if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
            lastwhite = whitechar(c)
        end do

        print *, "Number of columns in header = ", ncol

        rewind(unit)
        ! do
        !     read(unit,'(a)',iostat=ios) c
        !     if (iachar(c) == 37) then
        !         print *, "Skipped line"
        !         cycle
        !     else
        !         backspace(unit)
        !         exit
        !     end if
        ! end do    

        ! Parse values in header line
        ncon = 1
        cfmt = '000'
        select case(ncol)
        case(2)
            read(unit,*,iostat=ios) nvtxs, nedgs
        case(3)
            read(unit,*,iostat=ios) nvtxs, nedgs, cfmt
        case(4)
            read(unit,*,iostat=ios) nvtxs, nedgs, cfmt, ncon
            print *, "hello"
        case default
            write(*,*) "[load_graph]: incorrect file"
            stop
        end select
        read(cfmt,'(b3.3)') fmt
        print *, nvtxs, nedgs, cfmt, ncon
        write(*,'(A,B3.3)') "fmt = ", fmt

        ! Allocate necessary space
        allocate(xadj(nvtxs+1))
        allocate(adjncy(2*nedgs))

        if (btest(fmt,0)) allocate(adjwgt(2*nedgs))
        if (btest(fmt,1)) allocate(vwgt(nvtxs*ncon))
        if (btest(fmt,2)) allocate(vsize(nvtxs))

        write(*,*) associated(vsize),associated(vwgt),associated(adjwgt)
        ! stop


        xadj(1) = 0
        select case(fmt)
        case (b'000')
            do i = 1, nvtxs
                rowcol = count_columns(unit,stat=ios)
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) adjncy(xadj(i)+1:xadj(i+1))
            end do
        case(b'001')
            do i = 1, nvtxs
                rowcol = count_columns(unit,stat=ios)/2
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) (adjncy(j),adjwgt(j),j=xadj(i)+1,xadj(i+1))
            end do
        case(b'010')
            do i = 1, nvtxs
                rowcol = count_columns(unit,stat=ios) - ncon
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), adjncy(xadj(i)+1:xadj(i+1))
            end do
        case(b'100')
            do i = 1, nvtxs
                rowcol = count_columns(unit,stat=ios) - 1
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) vsize(i), adjncy(xadj(i)+1:xadj(i+1))
            end do
        case(b'011')
            do i = 1, nvtxs
                rowcol = (count_columns(unit,stat=ios) - ncon)/2
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), (adjncy(j),adjwgt(j),j=xadj(i)+1,xadj(i+1))
            end do
        case(b'110')
            do i = 1, nvtxs
                rowcol = count_columns(unit,stat=ios) - 1 - ncon
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) vsize(i), vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), adjncy(xadj(i)+1:xadj(i+1))
            end do
        case(b'101')
            do i = 1, nvtxs
                rowcol = (count_columns(unit,stat=ios) - 1)/2
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) vsize(i), (adjncy(j),adjwgt(j),j=xadj(i)+1,xadj(i+1))
            end do
        case(b'111')
            do i = 1, nvtxs
                rowcol = (count_columns(unit,stat=ios) - 1 - ncon)/2
                xadj(i+1) = xadj(i) + rowcol
                read(unit,*) vsize(i), vwgt((i-1)*ncon+1:(i-1)*ncon+ncon), (adjncy(j),adjwgt(j),j=xadj(i)+1,xadj(i+1))
            end do
        case default
            print *, "[read_graph] should not be here"
            stop
        end select

        if (numflag_ == 0) then
            adjncy = adjncy - 1
        else
            xadj = xadj + 1
        end if

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

    function FMETIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy,stat) result(ierr)
        use iso_c_binding, only : c_int, c_ptr, c_f_pointer
        integer(c_int), intent(in) :: ne
        integer(c_int), intent(in) :: nn
        integer(c_int), intent(in) :: eptr(ne+1)
        integer(c_int), intent(in) :: eind(:)
        integer(c_int), intent(in) :: numflag
        integer(c_int), intent(out), allocatable :: xadj(:)
        integer(c_int), intent(out), allocatable :: adjncy(:)
        integer(c_int), intent(out), optional :: stat ! stat = 0 indicates successful allocation and correct arguments

        ! Result
        integer(c_int) :: ierr

        integer(c_int) :: stat_
        character(len=80) :: errmsg_
        type(c_ptr) :: c_xadj, c_adjncy
        integer(c_int), pointer :: f_xadj(:) => null(), f_adjncy(:) => null()

        ierr = METIS_MeshToNodal(ne,nn,eptr,eind,numflag,c_xadj,c_adjncy)
        if (ierr /= METIS_OK) return

        call c_f_pointer(c_xadj,f_xadj,shape=[nn+1])
        
        select case(numflag)
        case(0)
            call c_f_pointer(c_adjncy,f_adjncy,shape=[f_xadj(nn+1)])
        case(1)
            call c_f_pointer(c_adjncy,f_adjncy,shape=[f_xadj(nn+1)-1])
        case default
            write(*,*) "[FMETIS_MeshToNodal] Wrong numflag argument! Only 0 or 1 are allowed. Got ", numflag, " instead."
            if (present(stat)) stat = -1
            return
        end select

        allocate(xadj,source=f_xadj,stat=stat_,errmsg=errmsg_)
        if (present(stat)) stat = stat_
        if (stat_ > 0) then
            write(*,*) "[FMETIS_MeshToNodal] Allocation of xadj failed with error: ", stat_, ", "//trim(errmsg_)//"."
            return
        end if

        allocate(adjncy,source=f_adjncy,stat=stat_,errmsg=errmsg_)
        if (present(stat)) stat = stat_
        if (stat_ > 0) then
            write(*,*) "[FMETIS_MeshToNodal] Allocation of adjncy failed with error: ", stat_, ", "//trim(errmsg_)//"."
            return
        end if

        ierr = METIS_Free(c_xadj)
        if (ierr /= METIS_OK) return

        ierr = METIS_Free(c_adjncy)
        if (ierr /= METIS_OK) return

    end function

end module
