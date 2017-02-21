#include "interval_partition.hpp"
#include "naive.hpp"
#include <gtest/gtest.h>

TEST(IntervalPartition, Example) {
	{
		using namespace IntervalPartition;
		IntervalPartition::IntervalledPolynom ip;
		{
			Polynom p(5);
			p[0] = 1;
			p[1] = Q(3,2);
			p[2] = Q(1,2);
			ip.push_back(9, std::move(p));
		}
		{
			Polynom p(5);
			p[0] = -44;
			p[1] = 11;
			ip.push_back(30, std::move(p));
		}
		{
			Polynom p(5);
			p[0] = -479;
			p[1] = Q(81,2);
			p[2] = Q(-1,2);
			ip.push_back(40, std::move(p));
		}
		{
			Polynom p(5);
			p[0] = 341;
			ip.push_back(49, std::move(p));
		}
		{
			Polynom p(5);
			p[0] = -884;
			p[1] = Q(99,2);
			p[2] = Q(-1,2);
			ip.push_back(59, std::move(p));
		}
		{
			Polynom p(5);
			p[0] = 946;
			p[1] = -11;
			ip.push_back(80, std::move(p));
		}
		{
			Polynom p(5);
			p[0] = 4186;
			p[1] = Q(-183,2);
			p[2] = Q(1,2);
			ip.push_back(90, std::move(p));
		}
		{
			const unsigned int dimensional_upper_bounds[] = {30, 50, 10};
			IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(dimensional_upper_bounds, sizeof(dimensional_upper_bounds)/sizeof(unsigned int), false);
			ASSERT_TRUE(intervalledPolynom == ip);
		}
		{
			unsigned int dimensional_upper_bounds[] = {10, 30, 50};
			const size_t dimensional_upper_bounds_length = sizeof(dimensional_upper_bounds)/sizeof(unsigned int);
			//TODO: accumulate
			const unsigned int maxdim = [&] () -> unsigned int
			{
				unsigned int maxdim_ = 0;
				for(const unsigned int& dim : dimensional_upper_bounds) maxdim_ += dim;
				return maxdim_;
			}();
			do {
				IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(dimensional_upper_bounds, dimensional_upper_bounds_length, false);
				for(size_t x = 0; x <= maxdim; ++x)
					ASSERT_EQ(intervalledPolynom(x), ip(x));
			}
			while(std::next_permutation(dimensional_upper_bounds, dimensional_upper_bounds+dimensional_upper_bounds_length));
		}
	}
}

#include "intpartrandom.hpp"

TEST_F(IntervalPartitionRandom, ModelCheck) {

	for(size_t steps = 0; steps < 1000; ++steps) {
		next();
		print();
		const auto naive_value = naive_bounds<mpz_class>(bounds, z, 0, bsize-1);
		const IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize, false);

		DVLOG(1) << naive_value;
		DVLOG(1) << intervalledPolynom.at(z);
		DVLOG(1) << intervalledPolynom(z);
		ASSERT_EQ(intervalledPolynom(z),  naive_value);
	}
}
TEST_F(IntervalPartitionRandom, PermutationInvariant) {
	for(size_t steps = 0; steps < 10; ++steps) {
		next();
		print();
		const IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize, false);
		const unsigned int maxdim = std::accumulate(bounds, bounds+bsize,0);

		// Invariant: If the bounds get shuffled, the resulting polynom's shape may change, but the resulting value should be the same.
		for(size_t permut_steps = 0; permut_steps < 100 && std::next_permutation(bounds, bounds+bsize); ++permut_steps) {
			IntervalPartition::IntervalledPolynom intervalledPolynom2 = IntervalPartition::generateIntervalPartition(bounds, bsize, false);
			for(size_t x = 0; x <= maxdim; ++x)
				ASSERT_EQ(intervalledPolynom(x), intervalledPolynom2(x));
		}
		// TODO: check if even the polynoms are equal!
	}
}

TEST_F(IntervalPartitionRandom, MirrorCheck) {
	for(size_t steps = 0; steps < 1000; ++steps) {
		next();
		print();
		IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize, false);
		IntervalPartition::IntervalledPolynom newPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize, true);
		const size_t maxdim = std::accumulate(bounds, bounds+bsize,static_cast<size_t>(0));
		DVLOG(1) << "O: " << intervalledPolynom.bounds().size() << "\n" << "N: " << newPolynom.bounds().size();

		for(size_t i = 0; i < maxdim/2; ++i)
			ASSERT_EQ(intervalledPolynom(i), newPolynom(i)) << intervalledPolynom << " vs " << newPolynom;
		//ASSERT_TRUE(intervalledPolynom == newPolynom) << intervalledPolynom << " vs " << newPolynom;
	}
}


TEST_F(IntervalPartitionRandom, ParallelCheck) {

	for(size_t steps = 0; steps < 1000; ++steps) {
		next();
		print();
		IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize,false);
		{
			IntervalPartition::IntervalledPolynom newPolynom = IntervalPartition::generateParallelIntervalPartition(bounds, bsize,false,1);
			ASSERT_TRUE(intervalledPolynom == newPolynom) << intervalledPolynom << " vs " << newPolynom;
		}
		{
			IntervalPartition::IntervalledPolynom newPolynom = IntervalPartition::generateParallelIntervalPartition(bounds, bsize,false,2);
			ASSERT_TRUE(intervalledPolynom == newPolynom) << intervalledPolynom << " vs " << newPolynom;
		}
		{
			IntervalPartition::IntervalledPolynom newPolynom = IntervalPartition::generateParallelIntervalPartition(bounds, bsize,false,3);
			ASSERT_TRUE(intervalledPolynom == newPolynom) << intervalledPolynom << " vs " << newPolynom;
		}
	}
}
