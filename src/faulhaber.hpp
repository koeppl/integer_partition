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
 *
 */
/**
 * @file faulhaber.hpp
 * @brief Pre-computes the Faulhaber polynoms
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
 */
#ifndef FAULHABER_HPP
#define FAULHABER_HPP
#include <cstddef>
#include "polynom.hpp"

namespace IntervalPartition
{

	/** Class for storing the Faulhaber polynoms
	 *
	 * The polynoms are generated in the constructor.
	*/
	class Faulhaber
	{
		public:
			const size_t dimension; //!< the number of polynoms that can be accessed by operator()
			
			/** 
			 * Precomputes the Faulhaber polynoms up to a certain dimension d
			 * @pre We need for the computation d+2 rows of Binomial numbers and
			 *      the first d+1 Bernoulli numbers
			 * 
			 * @param dimension The number of Faulhaber polynoms to compute. 
			 */
			Faulhaber(size_t dimension);
			~Faulhaber() {
				delete [] polynoms;
			}

			/** Accesses the i-th polynom
			 * 
			 * @param i \f$ 1 \le i \le \f$ dimension
			 * 
			 * @return the i-th Faulhaber polynom
			 */
			const Polynom& operator()(size_t i) const;
			const static Faulhaber f;
		private:
			Polynom*const polynoms;

	};

}//namespace
#endif//guard

