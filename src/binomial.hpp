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
 * @file binomial.hpp
 * @brief Pre-computes the binomial coefficients
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
 */
#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP
#include "definitions.hpp"

namespace IntervalPartition
{

	// TODO: exchange with mpz_bin_uiui (mpz_t rop, unsigned long int n, unsigned long int k) ?
	

	/** 
	 *  Stores the binomial coefficients in a two-dimensional array.
	 *  It is like Pascal's triangle, but the rows are "left-aligned".
	 */
	class Binomial
	{
		public:
			const size_t dimension; //<! the number of rows
			~Binomial() {
				for(size_t i = 0; i < dimension+1; ++i)
					delete [] binomials[i];
				delete binomials;
			}
			Binomial(size_t _dimension) 
				: dimension(_dimension), binomials(new Z*[dimension+1]), zero(0)
			{
				for(size_t i = 0; i < dimension+1; ++i) {
					binomials[i] = new Z[i+1];
					binomials[i][0] = 1;
				}
				for(size_t i = 1; i <= dimension; ++i)
				for(size_t j = 1; j <= i; ++j)
					binomials[i][j] = binomials[i-1][j-1] + binomials[i-1][j]; //TODO: binomials[][j] wrt. j symmetric! -> 1/2 space sufficient
			}
			const static Binomial b; //<! singleton instance

			/** 
			 * @param i must be less than dimension
			 * @param j any number is valid (out of bounds are catched by returning zero)
			 * 
			 * @return \$f i \choose j \$f
			 */
			const Z& operator()(size_t i, size_t j) const;
		private:
			Z**const binomials;
			const Z zero;
	};


}//namespace

std::ostream& operator<<(std::ostream& os, const IntervalPartition::Binomial& b);

#endif//guard
