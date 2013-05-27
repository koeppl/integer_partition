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
#ifndef POLYNOM_HPP
#define POLYNOM_HPP
#include "definitions.hpp"

namespace IntervalPartition
{

	class Polynom : public vektor<Q>
	{
		public:
		Polynom(size_t _size) : vektor<Q>(_size) {}
		Polynom(const Polynom& pol) : vektor<Q>(pol) {}
		Polynom() : vektor<Q>() {}
		Q operator()(const Q& x) const
		{
			Q res = at(size()-1);
			for(size_t i = 1; i < size(); ++i)
				res = res*x + at(size()-i-1);
			mpq_canonicalize(res.get_mpq_t());
			return res;
		}
		void canonicalize()
		{
			while(back() == 0) pop_back();
			for(Q& coeff : *this) mpq_canonicalize(coeff.get_mpq_t());
		}
		const static Polynom zero;
	};
	std::ostream& operator<<(std::ostream& os, const Polynom& v);
}
#endif//guard

