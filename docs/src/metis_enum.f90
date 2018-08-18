module metis_enum

    implicit none
    public

    !
    ! Enum type definitions
    !

    ! Return codes
    integer, parameter :: METIS_OK = 1              !! Returned normally.
    integer, parameter :: METIS_ERROR_INPUT = -2    !! Returned due to erroneous inputs and/or options.
    integer, parameter :: METIS_ERROR_MEMORY = -3   !! Returned due to insufficient memory.
    integer, parameter :: METIS_ERROR = -4          !! Return status: some other type of error.

    ! Operation type codes
    ! integer, parameter :: METIS_OP_PMETIS = 0       
    ! integer, parameter :: METIS_OP_KMETIS = 1
    ! integer, parameter :: METIS_OP_OMETIS = 2

    ! Options codes (i.e., `options`)
    integer, parameter :: METIS_OPTION_PTYPE     = 0  !! Specifies the partitioning method.
    integer, parameter :: METIS_OPTION_OBJTYPE   = 1  !! Specifies the type of objective.
    integer, parameter :: METIS_OPTION_CTYPE     = 2  !! Specifies the matching scheme to be used during coarsening.
    integer, parameter :: METIS_OPTION_IPTYPE    = 3  !! Determines the algorithm used during initial partitioning.
    integer, parameter :: METIS_OPTION_RTYPE     = 4  !! Determines the algorithm used for refinement.
    integer, parameter :: METIS_OPTION_DBGLVL    = 5  !! Specifies the amount of progress/debugging information will be printed.
    integer, parameter :: METIS_OPTION_NITER     = 6  !! Specifies the number of iterations for the refinement algorithm.
    integer, parameter :: METIS_OPTION_NCUTS     = 7  !! Specifies the number of different partitionings that it will compute.
    integer, parameter :: METIS_OPTION_SEED      = 8  !! Specifies the seed for the random number generator.
    integer, parameter :: METIS_OPTION_NO2HOP    = 9  !! Specifies that the coarsening will not perform any 2–hop matchings when the standard matching approach fails to sufficiently coarsen the graph.
    integer, parameter :: METIS_OPTION_MINCONN   = 10 !! Specifies that the partitioning routines should try to minimize the maximum degree of the subdomain graph.
    integer, parameter :: METIS_OPTION_CONTIG    = 11 !! Specifies that the partitioning routines should try to produce partitions that are contigous.
    integer, parameter :: METIS_OPTION_COMPRESS  = 12 !! Specifies that the graph should be compressed by combining together vertices that have identical adjacency lists.
    integer, parameter :: METIS_OPTION_CCORDER   = 13 !! Specifies if the connected components of the graph should first be identifies and ordered separately.
    integer, parameter :: METIS_OPTION_PFACTOR   = 14 !! Specifies the minimum degree of the vertices that will be ordered last.
    integer, parameter :: METIS_OPTION_NSEPS     = 15 !! Specifies the number of different separators that it will compute at each level of nested dissection.
    integer, parameter :: METIS_OPTION_UFACTOR   = 16 !! Specifies the maximum allowed load imbalance among the partitions.
    integer, parameter :: METIS_OPTION_NUMBERING = 17 !! Used to indicate which numbering scheme is used for the adjacency structure of a graph or the element-node structure of a mesh

    ! Partitioning Schemes
    integer, parameter :: METIS_PTYPE_RB   = 0 !! Multilevel recursive bisectioning.
    integer, parameter :: METIS_PTYPE_KWAY = 1 !! Multilevel \(k\)-way partitioning.

    ! Graph types for meshes
    ! integer, parameter :: METIS_GTYPE_DUAL  = 0
    ! integer, parameter :: METIS_GTYPE_NODAL = 1               

    ! Coarsening Schemes
    integer, parameter :: METIS_CTYPE_RM   = 0 !! Random matching.
    integer, parameter :: METIS_CTYPE_SHEM = 1 !! sorted heavy-edge matching.

    ! Initial partitioning schemes
    integer, parameter :: METIS_IPTYPE_GROW    = 0 !! Grows a bisection using a greedy strategy.
    integer, parameter :: METIS_IPTYPE_RANDOM  = 1 !! Computes a bisection at random followed by a refinement.
    integer, parameter :: METIS_IPTYPE_EDGE    = 2 !! Derives a separator form an edge cut.
    integer, parameter :: METIS_IPTYPE_NODE    = 3 !! Grows a bisection using a greedy node-based strategy.
    integer, parameter :: METIS_IPTYPE_METISRB = 4

    ! Refinement schemes
    integer, parameter :: METIS_RTYPE_FM        = 0 !! FM-based cut refinement.
    integer, parameter :: METIS_RTYPE_GREEDY    = 1 !! Greedy-based cut and volume refinement.
    integer, parameter :: METIS_RTYPE_SEP2SIDED = 2 !! Two-sided node FM refinement.
    integer, parameter :: METIS_RTYPE_SEP1SIDED = 3 !! One-sided node FM refinement.

    ! Debug Levels
    integer, parameter :: METIS_DBG_INFO       = 1       !! Shows various diagnostic messages.
    integer, parameter :: METIS_DBG_TIME       = 2       !! Perform timing analysis.
    integer, parameter :: METIS_DBG_COARSEN    = 4       !! Show the coarsening progress.
    integer, parameter :: METIS_DBG_REFINE     = 8       !! Show the refinement progress.
    integer, parameter :: METIS_DBG_IPART      = 16      !! Show info on initial partitioning.
    integer, parameter :: METIS_DBG_MOVEINFO   = 32      !! Show info on vertex moves during refinement.
    integer, parameter :: METIS_DBG_SEPINFO    = 64      !! Show info on vertex moves during sep refinement.
    integer, parameter :: METIS_DBG_CONNINFO   = 128     !! Show info on minimization of subdomain connectivity.
    integer, parameter :: METIS_DBG_CONTIGINFO = 256     !! Show info on elimination of connected components. 
    integer, parameter :: METIS_DBG_MEMORY     = 2048    !! Show info related to wspace allocation.
    
    ! Types of objectives
    integer, parameter :: METIS_OBJTYPE_CUT  = 0 !! Edge-cut minimization.
    integer, parameter :: METIS_OBJTYPE_VOL  = 1 !! Total communication volume minimization.
    integer, parameter :: METIS_OBJTYPE_NODE = 2

end module