# Source files

set( Numerics_SOURCES
	Numerics.h
	Numerics.cpp
)

include_directories(
	.
	../Common 
)


# Create the library
add_library(Numerics STATIC ${Numerics_SOURCES} )

#target_link_libraries( 
#	Numerics
#)

