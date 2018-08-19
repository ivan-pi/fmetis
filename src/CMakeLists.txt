
set(SRC
    metis_interface.f90
    metis_enum.f90
    metis_io.f90
    metis_oo_interface.f90
)

set(extern_INCLUDE /usr/local/include)
set(extern_LIB /usr/local/lib)
include_directories(${extern_INCLUDE})
link_directories (${extern_LIB})

add_library(fmetis ${SRC})
target_link_libraries(fmetis metis)
install(TARGETS fmetis
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
)