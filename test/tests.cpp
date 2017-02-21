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
#include "naive.hpp"


#include <random>
#include <gtest/gtest.h>
#include <gflags/gflags.h>


void singletest() {
	const unsigned int bounds[] = {2, 1, 2};
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int);
	constexpr size_t z = 0;

	std::cout << naive_bounds<mpz_class >(bounds, z, 0, bsize-1) << std::endl;
	IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize, false);
	std::cout << intervalledPolynom(z) << std::endl;
	DCHECK_EQ(intervalledPolynom(z),  naive_bounds<mpz_class >(bounds, z, 0, bsize-1));
}



// DEFINE_uint64(seed, std::default_random_engine::default_seed, "Random Seed");


namespace google {}
namespace gflags {}
int main(int argc, char **argv)
{
//	IB a = 25;
//	Q b = 4/7;
//	to_dump(a);
//	to_dump(b);
//	std::cout << to_string(a) << std::endl;
//	to_string(b);

	::testing::InitGoogleTest(&argc, argv);
	{
		using namespace google;
		using namespace gflags;
		ParseCommandLineFlags(&argc, &argv, true);
	}
	
	// if(FLAGS_seed != std::default_random_engine::default_seed)
	// 	RandomGenerator::m_seed = FLAGS_seed;
	// singletest();
	// return 0;
	return RUN_ALL_TESTS();
}


