#include "celero.hpp"
#include <gflags/gflags.h>
#include "binomial.hpp"
#include "faulhaber.hpp"
#include "bernoulli.hpp"
#include "interval_partition.hpp"
#include <celero/Celero.h>
#include "naive.hpp"

DEFINE_uint64(threads, 1, "Number of Threads");

int main(int argc, char **argv) {
	{
		using namespace google;
		using namespace gflags;
		ParseCommandLineFlags(&argc, &argv, true);
	}
		run_celero();
		return 0;
}

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
#define PAPERBENCH(s_number, s_bounds, s_z) \
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
	celero::DoNotOptimizeAway(IntervalPartition::generateParallelIntervalPartition(bounds, bsize, false, FLAGS_threads)(s_z)); \
}\
BENCHMARK(CONCATENATE(Roland, s_number), PartitionParallelMirror, 10, 10) \
{ \
	constexpr unsigned int bounds[] = s_bounds ; \
	constexpr size_t bsize = sizeof(bounds)/sizeof(unsigned int); \
	celero::DoNotOptimizeAway(IntervalPartition::generateParallelIntervalPartition(bounds, bsize, true, FLAGS_threads)(s_z)); \
} \















//PAPERBENCH( 12 , MACRO_ESCAPE({10000, 10005, 10010, 10015, 10021, 10027, 10039, 10063}), 40090)
//PAPERBENCH( 13 , MACRO_ESCAPE({10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000}), 5000)
//PAPERBENCH( 14 , MACRO_ESCAPE({10993, 10520, 10856, 10346, 10039, 10644, 10005, 10941, 10718, 10305}), 52683)
//PAPERBENCH( 15 , MACRO_ESCAPE({12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937}), 66131)
//PAPERBENCH( 19 , MACRO_ESCAPE({5641, 9314, 969, 8643, 6291, 6241, 8747, 7041}), 26433)
//PAPERBENCH( 17 , MACRO_ESCAPE({12184, 12324, 14685, 11098, 13357, 13863, 10796, 10914, 10989, 11115, 10937, 13634}), 72948)








PAPERBENCH( 1  , MACRO_ESCAPE({3, 4, 5}), 1)
PAPERBENCH( 2  , MACRO_ESCAPE({3, 4, 5}), 15)
PAPERBENCH( 3  , MACRO_ESCAPE({3000, 4000, 5000}), 1)
PAPERBENCH( 4  , MACRO_ESCAPE({3000, 4000, 5000}), 6000)
PAPERBENCH( 5  , MACRO_ESCAPE({10000, 10000 , 10000, 10000, 10000 }), 1)
PAPERBENCH( 6  , MACRO_ESCAPE({10000, 10000 , 10000, 10000, 10000 }), 300)
// PAPERBENCH( 7  , MACRO_ESCAPE({10000, 10000 , 10000, 10000, 10000}), 250000)
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


