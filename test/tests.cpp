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
	const unsigned int bounds[] = {2, 1, 2};
//	const unsigned int bounds[] = {18, 1, 79};
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int);
	constexpr size_t z = 0;

	std::cout << naive_bounds<mpz_class >(bounds, z, 0, bsize-1) << std::endl;
	IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize, false);
	std::cout << intervalledPolynom(z) << std::endl;
	DCHECK_EQ(intervalledPolynom(z),  naive_bounds<mpz_class >(bounds, z, 0, bsize-1));
}



#include "celero.hpp"
DEFINE_uint64(seed, std::default_random_engine::default_seed, "Random Seed");
DEFINE_bool(benchmark, false, "Do benchmark instead of testing");


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
	
//	if(FLAGS_seed != std::default_random_engine::default_seed)
//		RandomGenerator::m_seed = FLAGS_seed;
	if(FLAGS_benchmark) {
		run_celero();
		return 0;
	}
	//singletest();
	//return 0;
//	test_partition();
	return RUN_ALL_TESTS();
}


#include "binomial.hpp"
#include "faulhaber.hpp"
#include "bernoulli.hpp"
#include <celero/Celero.h>
BASELINE(Binomial, Baseline, 10, 10)
{
    celero::DoNotOptimizeAway(IntervalPartition::Binomial(BINOMIAL_DIM));
}
BASELINE(Bernoulli, Baseline, 10, 10)
{
    celero::DoNotOptimizeAway(IntervalPartition::Bernoulli(IntervalPartition::Binomial::b, BERNOULLI_DIM));
}
BASELINE(Faulhaber, Baseline, 10, 10)
{
    celero::DoNotOptimizeAway(IntervalPartition::Faulhaber(FAULHABER_DIM));
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

#include "macros.hpp"
#define ROLANDBENCH(s_number, s_bounds, s_z) \
BASELINE(CONCATENATE(Roland, s_number), Naive, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(naive_bounds<mpz_class >(bounds, s_z, 0, bsize-1)); \
} \
\
BENCHMARK(CONCATENATE(Roland, s_number), Partition, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateIntervalPartition(bounds, bsize, false)(s_z)); \
} \
\
BENCHMARK(CONCATENATE(Roland, s_number), PartitionMirror, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateIntervalPartition(bounds, bsize, true)(s_z)); \
} \
\
BENCHMARK(CONCATENATE(Roland, s_number), PartitionParallel, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateParallelIntervalPartition(bounds, bsize, false)(s_z)); \
}\
BENCHMARK(CONCATENATE(Roland, s_number), PartitionParallelMirror, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateParallelIntervalPartition(bounds, bsize, true)(s_z)); \
} 

















//ROLANDBENCH( 12 , MACRO_ESCAPE({10000, 10005, 10010, 10015, 10021, 10027, 10039, 10063}), 40090)
//ROLANDBENCH( 13 , MACRO_ESCAPE({10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000}), 5000)
//ROLANDBENCH( 14 , MACRO_ESCAPE({10993, 10520, 10856, 10346, 10039, 10644, 10005, 10941, 10718, 10305}), 52683)
//ROLANDBENCH( 15 , MACRO_ESCAPE({12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937}), 66131)
//ROLANDBENCH( 19 , MACRO_ESCAPE({5641, 9314, 969, 8643, 6291, 6241, 8747, 7041}), 26433)
//ROLANDBENCH( 17 , MACRO_ESCAPE({12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937, 13634}), 72948)








ROLANDBENCH( 1  , MACRO_ESCAPE({3, 4, 5}), 1)
ROLANDBENCH( 2  , MACRO_ESCAPE({3, 4, 5}), 15)
ROLANDBENCH( 3  , MACRO_ESCAPE({3000, 4000, 5000}), 1)
ROLANDBENCH( 4  , MACRO_ESCAPE({3000, 4000, 5000}), 6000)
ROLANDBENCH( 5  , MACRO_ESCAPE({104 , 104 , 104 , 104, 104 }), 1)
ROLANDBENCH( 6  , MACRO_ESCAPE({104 , 104 , 104, 104, 104 }), 300)
//ROLANDBENCH( 7  , MACRO_ESCAPE({104, 104, 104 , 104 , 104}), 250000)
//ROLANDBENCH( 8  , MACRO_ESCAPE({10000, 10005, 10010, 10015, 10020}), 25015)
ROLANDBENCH( 9  , MACRO_ESCAPE({10993, 10520, 10856, 10346, 10039}), 1)
//ROLANDBENCH( 10 , MACRO_ESCAPE({10993, 10520, 10856, 10346, 10039}), 26377)
ROLANDBENCH( 11 , MACRO_ESCAPE({33, 29, 42, 34, 59, 76, 54, 33}), 180)
//ROLANDBENCH( 16 , MACRO_ESCAPE({12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937, 13634}), 1)
ROLANDBENCH( 18 , MACRO_ESCAPE({3696, 3894, 4137, 7588, 7816}), 2856)
