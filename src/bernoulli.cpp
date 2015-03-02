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
#include "bernoulli.hpp"

namespace IntervalPartition
{
static_assert(BERNOULLI_DIM+1 < BINOMIAL_DIM, "Bernoulli dimension must be less than Binomial dimension -1");
const Bernoulli Bernoulli::b(Binomial::b, BERNOULLI_DIM);

std::ostream& operator<<(std::ostream& os, const IntervalPartition::Bernoulli& b) {
	os << "Bernoulli {";
	for(size_t i = 0; i < b.dimension; ++i) {
		os << b[i] << " ";
	}
	os << "}";
	os << std::endl;
	return os;
}

}//ns
