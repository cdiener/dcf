cmake_minimum_required(VERSION 2.8)

project(modes)

#Some Paths config
set(MODES_SRC ${CMAKE_SOURCE_DIR}/src)

include_directories(${MODES_SRC})
include_directories(${MODES_SRC}/alglib)
include_directories(${CMAKE_BINARY_DIR}/config)
add_subdirectory(${MODES_SRC}/alglib)

set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Curses)
if (CURSES_FOUND)
	include_directories(${CURSES_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${CURSES_LIBRARIES})
endif()

set(LIBS ${LIBS} alglib)

add_definitions(-std=c++11)

file(COPY examples DESTINATION .)
configure_file( ${MODES_SRC}/sys_config.h.in ${CMAKE_BINARY_DIR}/config/sys_config.h )

#Build program
add_executable(modes src/main.cxx src/iztli.cxx src/trainer.cxx src/design.cxx)
target_link_libraries(modes ${LIBS})

enable_testing()

#Test functionality
add_test (RunTest modes examples/alpha.txt 4 200 examples CPP efficiency)
set_tests_properties (RunTest
  PROPERTIES PASS_REGULAR_EXPRESSION "Needed.+for optimization\\.")

add_test(ArgFailTest modes bla 1)
set_tests_properties (ArgFailTest 
	PROPERTIES PASS_REGULAR_EXPRESSION "Usage: \\./modes")
