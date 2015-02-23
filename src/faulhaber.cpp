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
#include "faulhaber.hpp"
#include "binomial.hpp"
#include <cassert>
#include <glog/logging.h>

namespace IntervalPartition
{

/** 
 * Generates "naively" the Bernoulli numbers by using \class Binomial 
 * that pre-caches the Binomial coefficients
 */
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
		const static Bernoulli b;
	private:
		Q*const bernoulli;
};
const Bernoulli Bernoulli::b(Binomial::b, BERNOULLI_DIM);
}//namespace


std::ostream& operator<<(std::ostream& os, const IntervalPartition::Bernoulli& b) {
	os << "Bernoulli {";
	for(size_t i = 0; i < b.dimension; ++i) {
		os << b[i] << " ";
	}
	os << "}";
	os << std::endl;
	return os;
}

namespace IntervalPartition
{

const Faulhaber Faulhaber::f(FAULHABER_DIM);

Faulhaber::Faulhaber(size_t _dimension)
	: dimension(_dimension), polynoms(new Polynom[dimension])
{
	DCHECK_GT(BINOMIAL_DIM, dimension+2);
	DCHECK_GT(BERNOULLI_DIM, dimension+1);
	for(size_t p = 0; p < dimension; ++p)
	{
		polynoms[p].resize(p+2);
		for(size_t j = 0; j < p+1; ++j)
		{
			Q& coeff = polynoms[p][p+1-j];
			coeff = ( Bernoulli::b[j] * Binomial::b(p+1,j) )  / (p+1);
	//		std::cout << "(" << p+1 << "," << j << ") " << polynoms[j] << "\n";
			if(j%2) coeff *= -1;
		}
	}
}


const Polynom& Faulhaber::operator()(size_t i) const
{
	DCHECK_LT(i, dimension);
	return polynoms[i];
}

}//namespace
