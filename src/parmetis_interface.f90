
! http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download
! http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/manual.pdf

module parmetis_interface

    implicit none
    private

    ! Graph Partitioning
    public :: ParMETIS_V3_PartKway
    public :: ParMETIS_V3_PartGeomKway
    public :: ParMETIS_V3_PartGeom(
    public :: ParMETIS_V3_PartMeshKway

    ! Graph Repartitioning
    public :: ParMETIS_V3_AdaptiveRepart

    ! Partitioning Refinement
    public :: ParMETIS_V3_RefineKway

    ! Fill-reducing Orderings
    public :: ParMETIS_V3_NodeND
    public :: ParMETIS_V32_NodeND

    ! Mesh to Graph Translation
    public :: ParMETIS_V32_Mesh2Dual


    interface 

        !
        ! Graph Partitioning
        !

        function ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                      numflag,ncon,nparts,tpwgts,ubvec,options, &
                                      edgecut,part,comm) result(ierr) bind(C,name="ParMETIS_V3_PartKway")

        end function

        function ParMETIS_V3_PartGeomKway(vtxdist,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                          numflag,ndims,xyz,ncon,nparts,tpwgts, &
                                          ubvec,options,edgecut,part,comm) result(ierr) bind(C,name="ParMETIS_V3_PartGeomKway")

        end function

        function ParMETIS_V3_PartGeom(vtxdist,ndims,xyz,part,comm) result(ierr) bind(C,name="ParMETIS_V3_PartGeom")

        end function

        function ParMETIS_V3_PartMeshKway(elmdist,eptr,eind,elmwgt,wgtflag,numflag, &
                                          ncon,ncommon,nparts,tpwgts,ubvec, &
                                          options,edgecut,part,comm) result(ierr) bind(C,name="ParMETIS_V3_PartMeshKway")
        end function

        !
        ! Graph Repartitioning
        !

        function ParMETIS_V3_AdaptiveRepart(vtxdist,xadj,adjncy,vwgt,vsize,adjwgt, &
                                            wgtflag,numflag,ncon,nparts,tpwgts,ubvec, &
                                            itr,options,edgecut,part,comm) result(ierr) bind(C,name="ParMETIS_V3_AdaptiveRepart")

        !
        ! Partitioning Refinement
        !

        function ParMETIS_V3_RefineKway(vtxdist,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                        numflag,ncon,nparts,tpwgts,ubvec,options, &
                                        edgecut,part,comm) result(ierr) bind(C,name="ParMETIS_V3_RefineKway")

        !
        ! Fill-reducing Orderings
        !

        function ParMETIS_V3_NodeND(vtxdist,xadj,adjncy,numflag,options,order, &
                                    sizes,comm) result(ierr) bind(C,name="ParMETIS_V3_NodeND")
        end function

        function ParMETIS_V32_NodeND(vtxdist,xadj,adjncy,vwgt,numflag,mtype, &
                                     rtype,p_nseps,s_nseps,ubfrac,seed,dbglvl, &
                                     order,sizes,comm) result(ierr) bind(C,name="ParMETIS_V32_NodeND")

        end function

        !
        ! Mesh to Graph Translation
        !

        function ParMETIS_V32_Mesh2Dual(elmdist,eptr,eind,numflag,ncommon, &
                                        xadj,adjncy,comm) result(ierr) bind(C,name="ParMETIS_V32_Mesh2Dual")
            integer(c_int), intent(in) :: elmdist
            integer(c_int), intent(in) :: eptr, einf
            integer(c_int), intent(in) :: numflag
            integer(c_int), intent(in) :: ncommon
            type(c_ptr), intent(out) :: xadj,adjncy
            type(c_ptr), intent(in) :: comm
        end function

    end interface

end module