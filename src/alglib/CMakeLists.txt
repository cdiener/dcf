set(CMAKE_BUILD_TYPE Release)

add_library(alglib STATIC dataanalysis.h ap.h alglibinternal.h linalg.h statistics.h 
		specialfunctions.h alglibmisc.h solvers.h optimization.h 
		dataanalysis.cpp ap.cpp alglibinternal.cpp linalg.cpp statistics.cpp 
		specialfunctions.cpp alglibmisc.cpp solvers.cpp optimization.cpp)

