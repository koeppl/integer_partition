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
#include "polynom.hpp"

namespace IntervalPartition
{
	const Polynom Polynom::zero(1);
	const Polynom Polynom::one(1,1);
	std::ostream& operator<<(std::ostream& os, const Polynom& v)
	{
		os << "p{";
		for(size_t i = 0; i < v.size(); ++i)
		{
			os << v[i] << " ";
		}
		os << "}";
		return os;
	}
	Polynom operator-(const Polynom& a, const Polynom& b) {
		Polynom together(std::max(a.size(), b.size()));
		for(size_t l = 0; l < together.size(); ++l) {
			if(l < a.size()) together[l] += a[l];
			if(l < b.size()) together[l] -= b[l];
		}
		return together;
	}
	Polynom operator+(const Polynom& a, const Polynom& b) {
		Polynom together(std::max(a.size(), b.size()));
		for(size_t l = 0; l < together.size(); ++l) {
			if(l < a.size()) together[l] += a[l];
			if(l < b.size()) together[l] += b[l];
		}
		return together;
	}
}

