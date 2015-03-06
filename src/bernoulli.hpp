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
 * @file bernoulli.hpp
 * @brief Comuptation and storage of the Bernoulli Numbers
 * @author Dominik Köppl
 * 
 * @date 2015-03-02
 */

/** 
 * Generates "naively" the Bernoulli numbers by using \class Binomial 
 * that pre-caches the Binomial coefficients
 */
#include "binomial.hpp"
namespace IntervalPartition {

class Bernoulli
{
	public:
		const size_t dimension;
		Bernoulli(const Binomial& binomial, size_t _dimension)
			: dimension(_dimension), bernoulli(new Q[dimension])
		{
			DCHECK_GT(binomial.dimension, dimension+1);
			bernoulli[0] = Q(1,1);
			bernoulli[1] = Q(-1,2);
			for(size_t i = 2; i < dimension; ++i)
			{
				bernoulli[i] = 0;
				for(size_t j = 0; j < i; ++j) {
					bernoulli[i] -= bernoulli[j] * binomial(i+1,i+1-j);
				}
				bernoulli[i] /= i+1;
			}
		}
		~Bernoulli() {
			delete [] bernoulli;
		}
		Q operator[](size_t i) const {
			DCHECK_LT(i, dimension);
			return bernoulli[i];
		}
	const static Bernoulli b; //<! singleton instance
	private:
		Q*const bernoulli;
};
}//namespace

