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
/** 
 * @file sum_from_zero_to_upper.hpp
 * @brief Specialized Versions of the Faulhaber formula
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
 */

#ifndef SUM_FROM_ZERO_TO_UPPER
#define SUM_FROM_ZERO_TO_UPPER
#include <map>
#include "polynom.hpp"

namespace IntervalPartition
{   


	/**
	 * Computes the Faulhaber summation over a given polynom, cf. Lemma 4.5
	 */
	Polynom sumFromZeroToUpper(const Polynom& p);

/**
 * Calculates the polynom 
 * \f$ \sigma_{\gamma}(z) = \sum\limits_{k=0}^{z-\gamma} p(k) \f$
 * 
 * We use the function t_SumFunction to sum up the polynom 
 * and rewrite the coefficients
 * to match the \f$ z-\gamma \f$ upper bound of the sum.
 *
 * @tparam \function sumFromZeroToUpper or a function that caches these values
 * @tparam \class Binomial, Binomial::b or a function that returns binomial coefficients
 */
template<class t_SumFunction, class t_Binomial>
inline Polynom sumFromZeroToZMinusGamma(const Polynom& p, const Z& gamma, const t_SumFunction& sumFunction, const t_Binomial& binomial) {
	Polynom ret(p.size()+1);
	const Polynom& summedUp = sumFunction(p);
	for(size_t m = 0; m < p.size()+1; ++m) // m: index of leibniz binomial formula
	{
		Q& coeff = ret[m];
		Z gammapot = 1;
		for(size_t l = m; l < p.size()+1; ++l) {
			coeff += summedUp[l] * binomial(l,m) * gammapot;
			gammapot *= -gamma;
		}
	}
	return ret;
}

/** 
 * Computes \f$ \sum_{k = z - \gamma}^{upper} p(k) \f$
 * 
 * @return the polynom resulting by the summation
 *
 * We use the function t_SumFunction to sum up the polynom 
 * and rewrite the coefficients.
 *
 * The first term of the proof from Theorem 4.7
 * It is computed by the complete sum to the upper bound minus the sum from $0$ to $z-\gamma$.
 *
 * @tparam \function sumFromZeroToUpper or a function that caches these values
 */
template<class t_SumFunction, class t_Binomial>
inline Polynom sumFromZMinusGammaToUpper(const Polynom& p, const Z& gamma, const Z& upper, const t_SumFunction& sumFunction, const t_Binomial& binomial) {
	const Q highest = sumFunction(p)(upper);
	Polynom diff = sumFromZeroToZMinusGamma<t_SumFunction, t_Binomial>(p,gamma+1, sumFunction, binomial);
	for(Q& coeff : diff)
		coeff = -coeff;
	diff[0] += highest;
	return diff;
}


}//ns
#endif//guard




