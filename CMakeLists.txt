cmake_minimum_required(VERSION 3.5)

#set the project name
project(bielectronics_lib)

#include_folder
include_directories(${PROJECT_SOURCE_DIR})

#Create symbolic links to data files used during execution
ADD_CUSTOM_TARGET(link_target ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/src/FinalCubicSplineData.dat
		 	${CMAKE_CURRENT_BINARY_DIR}/FinalCubicSplineData.dat
		  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/src/CubicSplineDataF1num.dat
		 	${CMAKE_CURRENT_BINARY_DIR}/CubicSplineDataF1num.dat)


#add the executable

add_executable(bilib.x src/main.cpp src/MatElements.cpp src/McMurchie.cpp
	src/RysDupuisKing.cpp src/BoysFunction.cpp scr/readDataFiles.cpp 
	src/Coefficients.cpp src/Geometry.cpp src/Basis.cpp)
