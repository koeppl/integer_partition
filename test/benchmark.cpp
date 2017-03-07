#include <gflags/gflags.h>
#include "binomial.hpp"
#include "faulhaber.hpp"
#include "bernoulli.hpp"
#include "interval_partition.hpp"
#include <celero/Celero.h>
#include "naive.hpp"

namespace google {}
namespace gflags {}

//DEFINE_uint64(threads, std::thread::hardware_concurrency() == 0 ? 1 : std::thread::hardware_concurrency(), "Number of Threads");
DEFINE_uint64(threads, 4, "Number of Threads");

CELERO_MAIN

BASELINE(Binomial, Baseline, 100, 100)
{
    celero::DoNotOptimizeAway(IntervalPartition::Binomial(BINOMIAL_DIM));
}
BASELINE(Bernoulli, Baseline, 10, 100)
{
    celero::DoNotOptimizeAway(IntervalPartition::Bernoulli(IntervalPartition::Binomial::b, BERNOULLI_DIM));
}
BASELINE(Faulhaber, Baseline, 10, 100)
{
    celero::DoNotOptimizeAway(IntervalPartition::Faulhaber(FAULHABER_DIM));
}

#include "util.hpp"

namespace IntervalPartition {
	size_t computeValidityIntervals(const unsigned int* const dimensional_upper_bounds, 
			const size_t dimensions, 
			bool useSymmetry
			)
	{
		DVLOG(2) << "Interval Partitioning started";
		for(size_t i = 0; i < dimensions; ++i) {
			DCHECK_GT(dimensional_upper_bounds[i], 0) << "Every dimensional upper bound has to be > 0";
		}
		// maxdim is only used if useSymmetry to decide the cut-off
		const size_t maxdim = std::accumulate(dimensional_upper_bounds, dimensional_upper_bounds+dimensions,static_cast<size_t>(0)); 
		vektor<IB> intervalbounds;

		intervalbounds.push_back(dimensional_upper_bounds[0]);

		for(size_t k = 1; k < dimensions; ++k)
		{
			DVLOG(2) << "k: " << k;

			vektor<IB> tmp_intervalbounds; //! in this array the interval bounds of the next round (k+1) will be stored
			vektor<IB> help_intervalbounds;
			const unsigned int& dimensional_upper_bound = dimensional_upper_bounds[k];

			//if(dimensional_upper_bound > 1) 
			{ // we do now want a help_intervalbounds-value of 0
				help_intervalbounds.push_back(dimensional_upper_bound - 1);
			}


			for(const IB& old_intervalbound : intervalbounds)
				help_intervalbounds.push_back(old_intervalbound + dimensional_upper_bound);


			
			/**
			 * i : intervalbounds[] index
			 * j : help_intervalbounds[] index
			 * witness_right_index, witness_left_index : intervalbounds[] index
			 * witness_right, witness_left : help_intervalbounds[] value
			 *
			 */
			for(size_t i = 0, j = 0, witness_left_index = 0, witness_right_index = 0; j < help_intervalbounds.size();)
			{

				const IB& intervalbound = i < intervalbounds.size() ?  intervalbounds[i] : Z_zero;
				const IB& help_intervalbound = help_intervalbounds[j];
				
				const IB& witness_right = get_witness(witness_right_index, intervalbounds);
				const IB& witness_left = get_witness(witness_left_index, intervalbounds);
				if(useSymmetry && witness_left >= maxdim/2) break;
				DCHECK_LE(witness_left_index, intervalbounds.size()+1);
				DCHECK_LE(witness_right_index, intervalbounds.size()+1);

				DVLOG(2) << "k: " << k << ", i: " << i << ", j: " << j;
				if(i < intervalbounds.size()) {
					/**
					 * We are examining intervalbounds[i] and help_intervalbounds[j] and pop
					 * that value which is smaller (popping: increment either i or j)
					 * If both values are the same, we increment both i and j.
					 * min(intervalbounds[i], help_intervalbounds[j]) is appended to tmp_intervalbounds
					 */
					if(intervalbound < help_intervalbound) {
						tmp_intervalbounds.push_back(intervalbound);
						if(witness_right < intervalbound) {
							++witness_right_index;
						}
						if(intervalbound > witness_left + dimensional_upper_bound ) {
							++witness_left_index;
						}
						DVLOG(2) << "Fall 1, addiere " << intervalbound;
						++i;
					}
					else if(intervalbound > help_intervalbound)
					{
						tmp_intervalbounds.push_back(help_intervalbound);
						if(witness_right < help_intervalbound) {
							++witness_right_index;
						}
						if(help_intervalbound > witness_left + dimensional_upper_bound) {
							++witness_left_index;
						}
						DVLOG(2) << "Fall 2, addiere " << help_intervalbound;
						++j;
					}
					else
					{
						tmp_intervalbounds.push_back(intervalbound);
						if(intervalbound > witness_right) {
							++witness_right_index;
						}
						if(intervalbound > witness_left + dimensional_upper_bound) {
							++witness_left_index;
						}
						DVLOG(2) << "Fall 3, addiere " << intervalbound;
						++i;
						++j;
					}
				} else {
					/** We have already taken every element of intervalbounds[]. Because help_intervalbounds[] has some larger values, these have to 
					 *  be examined:
					 */
					tmp_intervalbounds.push_back(help_intervalbound);
					if(help_intervalbound > witness_right && witness_right_index <= intervalbounds.size()) {
						++witness_right_index;
					}
					if(help_intervalbound > witness_left + dimensional_upper_bound) {
						++witness_left_index;
					}
					DVLOG(2) << "Fall 4, addiere " << help_intervalbound;
					++j;
				}
				DVLOG(2) << "Witness: " << "[" << witness_left_index << ", " << witness_right_index << "]";
				DVLOG(2) << "tmp_intervalbounds: " << tmp_intervalbounds;
				DVLOG(2) << "intervalbounds: " << intervalbounds;

			}
			intervalbounds.swap(tmp_intervalbounds);
			DVLOG(2) << "_old_intervals: " << intervalbounds;
		}
		return intervalbounds.size();
	}
}



#include "macros.hpp"
#define PAPERBENCH(s_number, s_bounds, s_z) \
BASELINE(CONCATENATE(Paper, s_number), Naive, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(naive_bounds<mpz_class >(bounds, s_z, 0, bsize-1)); \
} \
\
BENCHMARK(CONCATENATE(Paper, s_number), Partition, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateIntervalPartition(bounds, bsize, false)(s_z)); \
} \
\
BENCHMARK(CONCATENATE(Paper, s_number), PartitionMirror, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateIntervalPartition(bounds, bsize, true)(s_z)); \
} \
\
BENCHMARK(CONCATENATE(Paper, s_number), PartitionParallel, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateParallelIntervalPartition(bounds, bsize, false, FLAGS_threads)(s_z)); \
}\
BENCHMARK(CONCATENATE(Paper, s_number), PartitionParallelMirror, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateParallelIntervalPartition(bounds, bsize, true, FLAGS_threads)(s_z)); \
} \
\
BENCHMARK(CONCATENATE(Paper, s_number), ValidityInterval, 100, 100) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::computeValidityIntervals(bounds, bsize, false)); \
}\












PAPERBENCH( 1  , MACRO_ESCAPE({3, 4, 5}), 1)
PAPERBENCH( 2  , MACRO_ESCAPE({3, 4, 5}), 15)
PAPERBENCH( 3  , MACRO_ESCAPE({3000, 4000, 5000}), 1)
PAPERBENCH( 4  , MACRO_ESCAPE({3000, 4000, 5000}), 6000)
PAPERBENCH( 5  , MACRO_ESCAPE({10000, 10000 , 10000, 10000, 10000 }), 1)
PAPERBENCH( 6  , MACRO_ESCAPE({10000, 10000 , 10000, 10000, 10000 }), 300)
// PAPERBENCH( 7  , MACRO_ESCAPE({10000, 10000 , 10000, 10000, 10000}), 25000)
// PAPERBENCH( 8  , MACRO_ESCAPE({10000, 10005, 10010, 10015, 10020}), 25015)
PAPERBENCH( 9  , MACRO_ESCAPE({10993, 10520, 10856, 10346, 10039}), 1)
// PAPERBENCH( 10 , MACRO_ESCAPE({10993, 10520, 10856, 10346, 10039}), 26377)
PAPERBENCH( 11 , MACRO_ESCAPE({33, 29, 42, 34, 59, 76, 54, 33}), 180)
// PAPERBENCH( 12 , MACRO_ESCAPE({10000,10005,10010,10015, 10021,10027,10039,10063}), 40090)
// PAPERBENCH( 13 , MACRO_ESCAPE({10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000}), 5000)
// PAPERBENCH( 14 , MACRO_ESCAPE({	10993, 10520, 10856, 10346, 10039, 10644, 10005, 10941, 10718, 10305              }),  52683)
// PAPERBENCH( 15 , MACRO_ESCAPE({	12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937       }),  66131)
// PAPERBENCH( 16 , MACRO_ESCAPE({	12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937, 13634}), 1)
// PAPERBENCH( 17 , MACRO_ESCAPE({	12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937, 13634}), 72948)
// PAPERBENCH( 18 , MACRO_ESCAPE({3696, 3894, 4137, 7588, 7816}), 2856)
// PAPERBENCH( 19 , MACRO_ESCAPE({	5641, 9314, 969, 8643, 6291, 6241, 8747,7041}),26433)


