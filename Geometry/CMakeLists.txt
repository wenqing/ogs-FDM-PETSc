# Source files

set( GEO_SOURCES
	GeoEntity.h
	Point.h
	Polyline.h
	Point.cpp
	Polyline.cpp
)

include_directories(
	.
	../Common
	../Numerics
)


# Create the library
add_library(Geometry STATIC ${GEO_SOURCES} )

#target_link_libraries( 
#	Geometry
#	Common
#)

