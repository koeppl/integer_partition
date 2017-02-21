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
#include "intervalled_polynom.hpp"
#include "polynom.hpp"
#include <glog/logging.h>
#include "util.hpp"

namespace IntervalPartition {

std::ostream& operator<<(std::ostream& os, const IntervalledPolynom& ip) {
	os << "IP: ";
	for(size_t i = 0; i < ip.intervalbounds.size(); ++i) {
		os << "{" << ip.intervalbounds[i] << " => " << ip.polynoms[i] << "}";
		if(i+1 < ip.intervalbounds.size()) os << ", ";
	}
	return os;
}

bool IntervalledPolynom::operator==(IntervalledPolynom& o) {
	for(size_t i = 0; i < intervalbounds.size(); ++i) {
		if(intervalbounds != o.intervalbounds) {
			DVLOG(2) << "intervalbounds not the same: " << intervalbounds[i] << " <-> " << o.intervalbounds[i];
			return false;
		}
	}
	for(size_t i = 0; i < intervalbounds.size(); ++i) {
		polynoms[i].canonicalize();
		o.polynoms[i].canonicalize();
		if(polynoms[i] != o.polynoms[i]) {
			DVLOG(2) << "Polynoms not the same: " << polynoms[i] << " <-> " << o.polynoms[i];
			return false;
		}
	}
	return true;
}
const Polynom& IntervalledPolynom::at(const IB& point) const
{
	if(point < 0) return Polynom::zero;
	DVLOG(2) <<  "Search " << point << " in " << intervalbounds;
	vektor<IB>::const_iterator it = lower_bound(intervalbounds.begin(), intervalbounds.end(), point);
	if(it != intervalbounds.end() && *it >= point) {
		const size_t distance = std::distance(intervalbounds.begin(),it);
		DVLOG(2) << "Found " << intervalbounds[distance] << "->" << polynoms[distance];
		return polynoms[distance];
	}
	DVLOG(2) << "Search failed !";
	return Polynom::zero;
}

Q IntervalledPolynom::operator()(const Z& x) const
{
	return at(x)(x);
}
/*
void IntervalledPolynom::push_back(const IB& intervalbound, const Polynom& polynom)
{
	intervalbounds.push_back(intervalbound);
	polynoms.push_back(polynom);
}
*/
void IntervalledPolynom::push_back(const IB& intervalbound, Polynom&& polynom)
{
	intervalbounds.push_back(intervalbound);
	polynoms.push_back(std::move(polynom));
}
void IntervalledPolynom::swap(IntervalledPolynom& o)
{
	intervalbounds.swap(o.intervalbounds);
	polynoms.swap(o.polynoms);
}


}

