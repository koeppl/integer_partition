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
 * @brief Caches polynom summations \f$ \sum\limits_{k=0}^{z} p(k) \f$ for any polynom p.
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
	 * A simple map is used to optimize lookups of the same value.
	 */
	class SumFromZeroToUpper
	{
		private:
		std::map<Polynom, Polynom> cache; //!< caches already made queries


		public:
		/** 
		 * Calculates the polynom 
		 * \f$ \sigma_{\gamma}(z) = \sum\limits_{k=0}^{z} p \f$
		 */
		const Polynom& operator()(const Polynom& p);
		static SumFromZeroToUpper s; //!< Singleton class
	};

	/** 
	 * Calculates the polynom 
	 * \f$ \sigma_{\gamma}(z) = \sum\limits_{k=0}^{z-\gamma} p(k) \f$
	 */
	Polynom sumFromZeroToZMinusGamma(const Polynom& p, const Z& gamma);

}
#endif//guard

