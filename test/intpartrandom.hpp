#include <gtest/gtest.h>
#include <glog/logging.h>
#include <random>

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
		std::stringstream ss;
		for(size_t i = 0; i < bsize; ++i) 
			ss << bounds[i] << ", ";
		DVLOG(1) << ss << z;
	}
};

