/* Integer Partition
 * Computes the number of possible ordered integer partitions with upper bounds
 * Copyright (C) 2013 Dominik Köppl
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


/*
#include <cassert>
#include <vector>
#include <map>
#include <iostream>
#include <cstdint>
#include <stdint.h>
#include <cmath>
#include <algorithm>

*/
#include "interval_partition.hpp"
#include "naive.hpp"


// output handling
/*
	template<class T>
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
	{
		os << "(";
		for(size_t i = 0; i < v.size(); ++i)
		{
			os << v[i] << " ";
		}
		os << ")";
		return os;
	}
	template<class U, class V>
	std::ostream& operator<<(std::ostream& os, const std::pair<U, V>& pair)
	{
		os << "|";
		os << pair.first << " -> " << pair.second;
		os << "|";
		return os;
	}

*/

#include <random>
#include <gtest/gtest.h>


void singletest() {
	const unsigned int bounds[] = {2, 1, 3, 1, 1, 1, 2, 1, 1};
//	const unsigned int bounds[] = {18, 1, 79};
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int);
	constexpr size_t z = 198;

	std::cout << naive_bounds<mpz_class >(bounds, z, 0, bsize-1) << std::endl;
	IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize);
	std::cout << intervalledPolynom(z) << std::endl;
	DCHECK_EQ(intervalledPolynom(z),  naive_bounds<mpz_class >(bounds, z, 0, bsize-1));
}

void test_partition() //(int argc, char** argv)
{
	std::default_random_engine generator;
	std::uniform_int_distribution<int> dim_distro(1,8); 
	std::uniform_int_distribution<int> urn_distro(2,50);

	for(int steps = 0; steps < 10000; ++steps) {
		const size_t bsize = dim_distro(generator);
		unsigned int*const bounds = new unsigned int[bsize];
		for(size_t i = 0; i < bsize; ++i)
			bounds[i] = urn_distro(generator);
		unsigned long z = std::uniform_int_distribution<unsigned long>(0, std::accumulate(bounds, bounds+bsize,0UL)+1)(generator);
		for(size_t i = 0; i < bsize; ++i) 
			std::cout << bounds[i] << ", ";
		std::cout << z << std::endl;
		std::cout << naive_bounds<mpz_class >(bounds, z, 0, bsize-1) << std::endl;
		IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize);
		std::cout << intervalledPolynom(z) << std::endl;
		EXPECT_EQ(intervalledPolynom(z),  naive_bounds<mpz_class >(bounds, z, 0, bsize-1));
		delete [] bounds;
	}


	/*
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
*/

}

#include "celero.hpp"
DEFINE_uint64(seed, std::default_random_engine::default_seed, "Random Seed");
DEFINE_bool(benchmark, false, "Do benchmark instead of testing");


namespace google {}
namespace gflags {}
int main(int argc, char **argv)
{
	IB a = 25;
	Q b = 4/7;
	to_dump(a);
	to_dump(b);
//	std::cout << to_string(a) << std::endl;
//	to_string(b);

	::testing::InitGoogleTest(&argc, argv);
	{
		using namespace google;
		using namespace gflags;
		ParseCommandLineFlags(&argc, &argv, true);
	}
	
//	if(FLAGS_seed != std::default_random_engine::default_seed)
//		RandomGenerator::m_seed = FLAGS_seed;
	if(FLAGS_benchmark) {
		run_celero();
		return 0;
	}
//	singletest();
//	return 0;
	test_partition();
	return 0;
	return RUN_ALL_TESTS();
}


#include "binomial.hpp"
#include <celero/Celero.h>
BASELINE(DemoSimple, Baseline, 10, 10)
{
    celero::DoNotOptimizeAway(IntervalPartition::Binomial(1000));
}
/*
BENCHMARK(DemoSimple, GMP, 10, 10)
{
   celero::DoNotOptimizeAway(IntervalPartition::Binomial(1000,true));
}
*/

//	long sum = 0;
//	for(size_t i = 0; i < bsize;++i) sum += bounds[i];
//	std::cout << intervalledPolynom(sum/2) << std::endl;




