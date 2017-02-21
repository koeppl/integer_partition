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


Q num_int_partition(unsigned int* const bounds, size_t bsize, unsigned long z) {
	size_t ones = 0; // number of dimensions with size 1
	for(size_t i = 0; i < bsize; ++i) {
		if(bounds[i] == 0) return 0;
		if(bounds[i] == 1) {
			bounds[i] = bounds[bsize-1];
			++ones;
			--bsize;
			--i;
		}
	}
	const size_t dimensionalSum = std::accumulate(bounds, bounds+bsize, static_cast<size_t>(0));
	IntervalPartition::IntervalledPolynom intervalledPolynom = FLAGS_threads == 1
		? IntervalPartition::generateIntervalPartition(bounds, bsize, true)
		: IntervalPartition::generateParallelIntervalPartition(bounds, bsize, true, FLAGS_threads);
	if(z > dimensionalSum) {
		z  = dimensionalSum-z;
	}
	if(ones == 0) { 
		return intervalledPolynom(z);
	}
	Q ret = 0;
	const size_t sum_bound = std::min(z, ones);
	for(size_t k = 0; k <= sum_bound; ++k) {
		ret += IntervalPartition::Binomial::b(ones,k) * intervalledPolynom(z-k);
	}
	return ret;
}

int main(int argc, char** argv) {
	{
		using namespace google;
		using namespace gflags;
		ParseCommandLineFlags(&argc, &argv, true);
	}
	if(argc < 3)
	{
		std::cout << argv[0] << " - calculate the " << std::endl;
		std::cout << "Usage: " << argv[0] << " n i_0 [i_1 [i_2 [...]]]" << std::endl;
		return 1;
	}
	const size_t bsize = argc-2;
	const unsigned long z = strtoul(argv[1], NULL, 10);
    unsigned int*const bounds = new unsigned int[bsize];
	for(size_t i = 2; i < static_cast<size_t>(argc); ++i)
		bounds[i-2] = strtoul(argv[i], NULL, 10);

	std::cout << num_int_partition(bounds, bsize, z) << std::endl;

	delete [] bounds;
	return 0;
}

//	long sum = 0;
//	for(size_t i = 0; i < bsize;++i) sum += bounds[i];
//	std::cout << intervalledPolynom(sum/2) << std::endl;



