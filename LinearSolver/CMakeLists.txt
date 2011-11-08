# Source files


IF(USE_PETSC)
	set( LinSolver_SOURCES
	     PETScLinearSolver.h
	     PETScLinearSolver.cpp
	)
ELSE(USE_PETSC)
	set( LinSolver_SOURCES
	     LinearEQS.h
	     LinearEQS.cpp
	)
ENDIF(USE_PETSC)

include_directories(
	.
	../Common
	../Math
	../Numerics
)


# Create the library
add_library(LinearSolver STATIC ${LinSolver_SOURCES} )


IF(USE_PETSC)
	target_link_libraries( 
	  LinearSolver
	)
ELSE(USE_PETSC)
	target_link_libraries( 
	  LinearSolver
	  Math
	)
ENDIF(USE_PETSC)
