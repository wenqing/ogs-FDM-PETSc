# Source files

set( Common_SOURCES
	misc.h
	misc.cpp
)

include_directories(
	.
)


# Create the library
add_library(Common STATIC ${Common_SOURCES} )

#target_link_libraries( 
#	Common
#)

