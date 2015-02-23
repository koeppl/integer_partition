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

include(CMakeListsFiles.txt)
add_library(integer_partition SHARED ${integer_partition_SRCS})
target_link_libraries(integer_partition ${GMP_LIBRARY} ${GMPXX_LIBRARY})
target_link_libraries(integer_partition ${GLOG_LIBRARY})

add_library(integer_partition_static STATIC ${integer_partition_SRCS})
set_target_properties(integer_partition_static PROPERTIES OUTPUT_NAME integer_partition)

INSTALL_TARGETS(/lib integer_partition)
INSTALL_TARGETS(/lib integer_partition_static)
INSTALL_FILES(/include/integer_partition .hpp)
