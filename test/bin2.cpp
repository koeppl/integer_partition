#include <celero/Celero.h>
#include "binomial.hpp"
#include <gtest/gtest.h>

	class Binomial2
	{
		public:
			const size_t dimension; //<! the number of rows
			~Binomial2() {
				for(size_t i = 0; i < dimension+1; ++i)
					delete [] binomials[i];
				delete binomials;
			}
			Binomial2(size_t _dimension) 
				: dimension(_dimension), binomials(new Z*[dimension+1]) 
			{
				for(size_t i = 0; i <= dimension; ++i) {
					binomials[i] = new Z[i/2+2]; //We need i+2 as we use later binomials[i-1][j]
					binomials[i][0] = 1;
					//binomials[i][i+1] = 0; //!< reset to zero if the default constructor of Z does not do this job
				}
				for(size_t i = 1; i <= dimension; ++i) {
					for(size_t j = 1; j <= i/2; ++j) {
						binomials[i][j] = binomials[i-1][j-1] + binomials[i-1][j]; //TODO: binomials[][j] wrt. j symmetric! -> 1/2 space sufficient
//						std::cout << i << "," << j << "->" << binomials[i][j] << std::endl;
					}
					if(i % 2) {
					binomials[i][i/2+1] = 2*binomials[i-1][(i-1)/2];
//						std::cout << i << "," << i/2+1 << "->" << binomials[i][i/2+1] << std::endl;
					}
				}
			
			}

			/** 
			 * @param i must be less than dimension
			 * @param j any number is valid (out of bounds are catched by returning zero)
			 * 
			 * @return \$f i \choose j \$f
			 */
			const Z& operator()(size_t i, size_t j) const;
		private:
			Z**const binomials;
	};
const Z& Binomial2::operator()(size_t i, size_t j) const {
	DCHECK_LT(i,  dimension);
	if(i == 0) return Z_zero;
	if(j > i) return Z_zero;
	if(j > i/2) return binomials[i-1][i-j];
	return binomials[i-1][j];
}

BENCHMARK(Binomial, BinTwo, 10, 10)
{
    celero::DoNotOptimizeAway(Binomial2(BINOMIAL_DIM));
}
BASELINE(BinomialAccess, Baseline, 10, 100)
{
	for(size_t i = 0; i < IntervalPartition::Binomial::b.dimension; ++i)
		for(size_t j = 0; j < i; ++j)
			celero::DoNotOptimizeAway(IntervalPartition::Binomial::b(i,j));
}
BENCHMARK(BinomialAccess, BinTwo, 10, 100)
{
	Binomial2 bin(BINOMIAL_DIM);
	for(size_t i = 0; i < IntervalPartition::Binomial::b.dimension; ++i)
		for(size_t j = 0; j < i; ++j)
			celero::DoNotOptimizeAway(bin(i,j));
}

TEST(Binomial, BinTwo) {
	Binomial2 bin(BINOMIAL_DIM);
	for(size_t i = 1; i < IntervalPartition::Binomial::b.dimension; ++i)
//	Binomial2 bin(10);
//	for(size_t i = 0; i < 10; ++i)
	for(size_t j = 0; j < i; ++j) {
		Q c = bin(i,j);
		Q b = IntervalPartition::Binomial::b(i,j);
		mpq_canonicalize(c.get_mpq_t());
		ASSERT_TRUE(mpq_equal(c.get_mpq_t(), b.get_mpq_t())) << "(" << i << " choose " << j << "), " << c << " vs " << b;
	}	
}
