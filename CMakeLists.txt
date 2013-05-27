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
set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=gnu++11 -pedantic -Wshadow")


find_library(GMP_LIBRARY gmp /usr/lib)
find_library(GMPXX_LIBRARY gmpxx /usr/lib)

add_subdirectory (src)
add_subdirectory (demo) 


ENABLE_TESTING()
include_directories (${integer_partition_SOURCE_DIR}/src) 
link_directories (${integer_partition_BINARY_DIR}/src) 
add_executable (integer_partition_test test.cpp)
target_link_libraries (integer_partition_test integer_partition)
ADD_TEST(integer_partition_test ${CMAKE_CURRENT_BINARY_DIR}/integer_partition_test)
