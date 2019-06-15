include_directories(${PROJECT_BINARY_DIR}/include)
link_directories(${CMAKE_BINARY_DIR}/lib)

add_executable(oo_tests oo_tests.f90)
target_link_libraries(oo_tests fmetis)