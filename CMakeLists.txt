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

#my own profile
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /home/niki/opt/lib)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /home/niki/opt/include)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I/home/niki/opt/include" )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L/home/niki/opt/lib" )



set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

project(integer_partition)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=gnu++11 -pedantic -O0 -ggdb -Wshadow")


find_library(GMP_LIBRARY gmp /usr/lib)
find_library(GMPXX_LIBRARY gmpxx /usr/lib)
find_package(Glog)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})
find_package (GFlags REQUIRED)
include_directories (${gflags_INCLUDE_DIR})


add_subdirectory (test)
add_subdirectory (src)
add_subdirectory (demo) 
add_subdirectory (naive)


ENABLE_TESTING()
include_directories (${integer_partition_SOURCE_DIR}/src) 
link_directories (${integer_partition_BINARY_DIR}/src) 
add_executable (integer_partition_test test.cpp)
target_link_libraries (integer_partition_test integer_partition)
ADD_TEST(integer_partition_test ${CMAKE_CURRENT_BINARY_DIR}/integer_partition_test)
