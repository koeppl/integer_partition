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
#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP
#include "definitions.hpp"

namespace IntervalPartition
{

	class Binomial
	{
		public:
			const size_t dimension;
			~Binomial()
			{
				for(size_t i = 0; i < dimension+1; ++i)
					delete [] binomials[i];
				delete binomials;
			}
			Binomial(size_t _dimension) 
				: dimension(_dimension), binomials(new Z*[dimension+1]), zero(0)
			{
				for(size_t i = 0; i < dimension+1; ++i)
				{
					binomials[i] = new Z[i+1];
					binomials[i][0] = 1;
				}
				for(size_t i = 1; i <= dimension; ++i)
				for(size_t j = 1; j <= i; ++j)
					binomials[i][j] = binomials[i-1][j-1] + binomials[i-1][j]; //TODO: binomials[][j] wrt. j symmetric! -> 1/2 space sufficient
			}
			const static Binomial b;
			const Z& operator()(size_t i, size_t j) const;
		private:
			Z**const binomials;
			const Z zero;
	};


}//namespace

std::ostream& operator<<(std::ostream& os, const IntervalPartition::Binomial& b);

#endif//guard
