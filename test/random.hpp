/** 
 * @file random.hpp
 * @brief Wrapper around libGMP for creating random variables.
 * @author Dominik Köppl
 * 
 * @date 2015-03-02
 */

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "gmpxx.h"
namespace IntervalPartition {

/** 
 * Wrapper around libGMP for creating random variables.
 * Returns random variables of type \class mpz_class
 */
class GMPRandom
{
	gmp_randstate_t r_state;
	public:

	/** 
	 * @param seed the seed with which to initialize
	 */
	GMPRandom(unsigned long seed) {
		gmp_randinit_default(r_state);
		gmp_randseed_ui(r_state, seed);
	}
	~GMPRandom() {
		gmp_randclear(r_state);
	}
	mpz_class get() {
		mpz_class z;
		mpz_urandomb(z.get_mpz_t(),r_state,8);
		return z;
	}
};

}//ns
#endif /* RANDOM_HPP */
