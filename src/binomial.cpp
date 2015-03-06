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

#include "binomial.hpp"
#include "glog/logging.h"

namespace IntervalPartition {



const Z& Binomial::operator()(size_t i, size_t j) const {
	DCHECK_LT(i,  dimension);
	if(j > i) return Z_zero;
	return binomials[i][j];
}

}//namespace

std::ostream& operator<<(std::ostream& os, const IntervalPartition::Binomial& b) {
	for(size_t i = 0; i < b.dimension; ++i) {
		os << "(";
		for(size_t j = 0; j <= i; ++j)
			os << b(i,j) << " ";
		os << ")\n";
	}
	os << std::endl;
	return os;
}
