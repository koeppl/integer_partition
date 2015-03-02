# Integer Partition
# Computes the number of possible ordered integer partitions with upper bounds
# Copyright (C) 2013 Dominik Köppl
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.


cmake_minimum_required(VERSION 2.8)
project(integer_partition)

#my own profile
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /home/niki/opt/lib)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /home/niki/opt/include)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I/home/niki/opt/include" )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L/home/niki/opt/lib" )


set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=gnu++11 -pedantic -ggdb -Wshadow")


find_library(GMP_LIBRARY gmp /usr/lib)
find_library(GMPXX_LIBRARY gmpxx /usr/lib)
find_package(Glog)
find_package(Threads)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})
find_package (GFlags REQUIRED)
include_directories (${gflags_INCLUDE_DIR})



option(CODE_COVERAGE "Set ON to add code coverage compile options" OFF)
if( CODE_COVERAGE )
	set(CMAKE_BUILD_TYPE "Debug")
	include(CodeCoverage)
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage" )
	SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lgcov" )
endif()

option(CODE_PROFILING "Set ON to add profiling compile options" OFF)
if( CODE_PROFILING )
	set(CMAKE_BUILD_TYPE "Release")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pg" )
endif()

if (CMAKE_BUILD_TYPE EQUAL "Release")
	add_definitions(-DGOOGLE_STRIP_LOG=1)
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

MESSAGE( STATUS "With Code Coverage?: " ${CODE_COVERAGE} )
MESSAGE( STATUS "With Code Profiling?: " ${CODE_PROFILING} )

MESSAGE( STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CXX Flags: " ${CMAKE_CXX_FLAGS} )




add_subdirectory (src)
add_subdirectory (demo) 
add_subdirectory (test)
add_subdirectory (naive)


#	ENABLE_TESTING()
#	include_directories (${integer_partition_SOURCE_DIR}/src) 
#	link_directories (${integer_partition_BINARY_DIR}/src) 
#	add_executable (integer_partition_test test.cpp)
#	target_link_libraries (integer_partition_test integer_partition)
#	ADD_TEST(integer_partition_test ${CMAKE_CURRENT_BINARY_DIR}/integer_partition_test)
