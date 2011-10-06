# Source files

set(FDM_SOURCES
ConditionDataBC.cpp
ConditionDataBC.h
FiniteDifference.cpp
FiniteDifference.h
MatProperty.cpp
MatProperty.h
Output.cpp
Output.h
RasterRecharge.cpp
RasterRecharge.h
)

include_directories(
	.
	../Math 
	../Numerics 
	../Geometry 
    ../LinearSolver  
    ../Common   
    ../Numerics   
)


# Create the library


add_library(FDM STATIC ${FDM_SOURCES} )

add_executable(${PROJECT_NAME} main.cpp)


target_link_libraries(${PROJECT_NAME} 
FDM
Math
Geometry	
LinearSolver
Common
Numerics
)

