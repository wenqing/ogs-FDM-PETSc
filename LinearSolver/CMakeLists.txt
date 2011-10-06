# Source files

set( LinSolver_SOURCES
	LinearEQS.h
	LinearEQS.cpp
)

include_directories(
	.
	../Math
	../Numerics
)


# Create the library
add_library(LinearSolver STATIC ${LinSolver_SOURCES} )

#target_link_libraries( 
#	LinearSolver
#)

