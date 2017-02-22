/* Integer Partition
 * Computes the number of possible ordered integer partitions with upper bounds
 * Copyright (C) 2013,2017 Dominik Köppl
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "interval_partition.hpp"
#include "binomial.hpp"
#include <gflags/gflags.h>

DEFINE_uint64(threads, 1, "Number of Threads");


namespace gflags {}
namespace google {}

int main(int argc, char** argv) {
	{
		using namespace google;
		using namespace gflags;
		ParseCommandLineFlags(&argc, &argv, true);
	}
	if(argc < 3)
	{
		std::cout << argv[0] << " - calculate the " << std::endl;
		std::cout << "Usage: " << argv[0] << " z i_0 [i_1 [i_2 [...]]]" << std::endl;
		return 1;
	}
	const size_t bsize = argc-2;
	const unsigned long z = strtoul(argv[1], NULL, 10);
    unsigned int*const bounds = new unsigned int[bsize];
	for(size_t i = 2; i < static_cast<size_t>(argc); ++i)
		bounds[i-2] = strtoul(argv[i], NULL, 10);

	std::cout << IntervalPartition::number_of_interval_partitions(bounds, bsize, z, FLAGS_threads) << std::endl;

	delete [] bounds;
	return 0;
}



