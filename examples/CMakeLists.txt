include_directories(${PROJECT_BINARY_DIR}/src)

add_executable(metis_tests metis_tests.f90)
target_link_libraries(metis_tests fmetis)

add_executable(oo_tests oo_tests.f90)
target_link_libraries(oo_tests fmetis)