cmake_minimum_required(VERSION 2.8)

project(cppdesigner)

#Some Paths config
set(CPPDESIGNER_SRC ${CMAKE_SOURCE_DIR}/src)

include_directories(${CPPDESIGNER_SRC})
include_directories(${CPPDESIGNER_SRC}/alglib)
add_subdirectory(${CPPDESIGNER_SRC}/alglib)

set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_definitions(-std=c++0x)

file(COPY examples DESTINATION .)

#Build program
add_executable(cppdesigner src/main.cxx src/iztli.cxx src/design.cxx)
target_link_libraries(cppdesigner alglib)

enable_testing()

#Test functionality
add_test (RunTest cppdesigner examples/alpha.txt 4 200 examples/pos.txt examples/neg.txt examples/high.txt examples/low.txt)
set_tests_properties (RunTest
  PROPERTIES PASS_REGULAR_EXPRESSION "error estimate.+Link 0.+optimization\\.")

add_test(ArgFailTest cppdesigner bla 1)
set_tests_properties (ArgFailTest 
	PROPERTIES PASS_REGULAR_EXPRESSION "Usage: ./cppdesigner sequence_file n_linker n_iter")
