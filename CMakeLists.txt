cmake_minimum_required(VERSION 2.8)

project(dcf)

#Some Paths config
set(DCF_SRC ${CMAKE_SOURCE_DIR}/src)

include_directories(${DCF_SRC})
include_directories(${DCF_SRC}/alglib)
include_directories(${CMAKE_BINARY_DIR}/config)
add_subdirectory(${DCF_SRC}/alglib)

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

set(LIBS ${LIBS} alglib -Wl,-O1,--sort-common,--as-needed,-z,relro)

add_definitions(-std=c++11 -pipe -fstack-protector --param=ssp-buffer-size=4 -D_FORTIFY_SOURCE=2)

file(COPY examples DESTINATION .)
file(COPY scripts DESTINATION .)
file(COPY src/blosum80.qij DESTINATION .)
configure_file( ${DCF_SRC}/sys_config.h.in ${CMAKE_BINARY_DIR}/config/sys_config.h )

#Build program
add_executable(dcf src/main.cxx src/iztli.cxx src/trainer.cxx src/design.cxx)
target_link_libraries(dcf ${LIBS})

#Build predict
add_executable(predict src/predict.cxx src/iztli.cxx src/trainer.cxx)
target_link_libraries(predict ${LIBS})

#Build pfam_splitter
add_executable(pfam_splitter src/split.cxx src/pfam.cxx)
target_link_libraries(pfam_splitter ${LIBS})

#Build pfam_splitter
add_executable(prop src/prop.cxx src/iztli.cxx)
target_link_libraries(pfam_splitter ${LIBS})

#Build tests
add_executable(tests src/tests.cxx src/iztli.cxx src/trainer.cxx src/design.cxx)
target_link_libraries(tests ${LIBS})
