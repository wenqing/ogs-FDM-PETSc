# Source files

set( Math_SOURCES
    Matrix.h
    Matrix.cpp    
	SparseTable.h
	SparseTable.cpp
	SparseMatrix.h
	SparseMatrix.cpp
)

include_directories(
	.
	../Math
	../FDM
	../Common	
	../Geometry
)


# Create the library
add_library(Math STATIC ${Math_SOURCES} )

#target_link_libraries( 
#	Math
#)

