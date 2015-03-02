#include "interval_partition.hpp"
#include "naive.hpp"
#include <gtest/gtest.h>
#include <random>

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
			IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(dimensional_upper_bounds, sizeof(dimensional_upper_bounds)/sizeof(unsigned int));
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
				IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(dimensional_upper_bounds, dimensional_upper_bounds_length);
				for(size_t x = 0; x <= maxdim; ++x)
					ASSERT_EQ(intervalledPolynom(x), ip(x));
			}
			while(std::next_permutation(dimensional_upper_bounds, dimensional_upper_bounds+dimensional_upper_bounds_length));
		}
	}
}

class IntervalPartitionRandom : public ::testing::Test {
	std::default_random_engine generator;
	std::uniform_int_distribution<int> dim_distro = std::uniform_int_distribution<int>(1,10); 
	std::uniform_int_distribution<int> urn_distro = std::uniform_int_distribution<int>(2,10);

	protected:
	unsigned int* bounds = nullptr;
	size_t bsize;
	unsigned long z;

	virtual void SetUp() override {
		bounds = nullptr;
	}
	virtual void TearDown() {
		if(bounds != nullptr) delete [] bounds;
	}


	void next() {
		bsize = dim_distro(generator);
		if(bounds != nullptr) delete [] bounds;
		bounds = new unsigned int[bsize];
		for(size_t i = 0; i < bsize; ++i)
			bounds[i] = urn_distro(generator);
		z = std::uniform_int_distribution<unsigned long>(0, std::accumulate(bounds, bounds+bsize,0UL)+1)(generator);
	}
	unsigned long get_z() const { return z; }
	unsigned int* get_bounds() { return bounds; }
	size_t get_bsize() const { return bsize; }

	void print() const {
		for(size_t i = 0; i < bsize; ++i) 
			std::cout << bounds[i] << ", ";
		std::cout << z << std::endl;
	}
};

TEST_F(IntervalPartitionRandom, ModelCheck) {

	for(size_t steps = 0; steps < 1000; ++steps) {
		next();
		print();
		const auto naive_value = naive_bounds<mpz_class>(bounds, z, 0, bsize-1);
		const IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize);

		std::cout << naive_value << std::endl;
		std::cout << intervalledPolynom.at(z) << std::endl;
		std::cout << intervalledPolynom(z) << std::endl;
		ASSERT_EQ(intervalledPolynom(z),  naive_value);
	}
}
TEST_F(IntervalPartitionRandom, PermutationInvariant) {
	for(size_t steps = 0; steps < 10; ++steps) {
		next();
		print();
		const IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize);
		const unsigned int maxdim = std::accumulate(bounds, bounds+bsize,0);

		// Invariant: If the bounds get shuffled, the resulting polynom's shape may change, but the resulting value should be the same.
		for(size_t permut_steps = 0; permut_steps < 100 && std::next_permutation(bounds, bounds+bsize); ++permut_steps) {
			IntervalPartition::IntervalledPolynom intervalledPolynom2 = IntervalPartition::generateIntervalPartition(bounds, bsize);
			for(size_t x = 0; x <= maxdim; ++x)
				ASSERT_EQ(intervalledPolynom(x), intervalledPolynom2(x));
		}
		// TODO: check if even the polynoms are equal!
	}
}
