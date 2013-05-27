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
#ifndef INTERVALLED_POLYNOM
#define INTERVALLED_POLYNOM
#include "definitions.hpp"
#include "polynom.hpp"

namespace IntervalPartition
{

	class Polynom;
	class IntervalledPolynom;
	std::ostream& operator<<(std::ostream& os, const IntervalledPolynom& ip);

	class IntervalledPolynom
	{
		private:
			vektor<IB> intervalbounds;
			vektor<Polynom> polynoms;
		public:
			const Polynom& at(const IB& point) const;
			void push_back(const IB& intervalbound, const Polynom& polynom);
			void swap(IntervalledPolynom& o);
			Q operator()(const Z& x) const;
			bool operator==(IntervalledPolynom& o);

		friend std::ostream& operator<<(std::ostream& os, const IntervalledPolynom& ip);
	};

}//namespace
std::ostream& operator<<(std::ostream& os, const IntervalPartition::IntervalledPolynom& ip);
#endif//guard

