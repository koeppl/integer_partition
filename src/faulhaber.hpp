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
#ifndef FAULHABER_HPP
#define FAULHABER_HPP
#include <cstddef>
#include "polynom.hpp"

namespace IntervalPartition
{

	class Faulhaber
	{
		public:
			const size_t dimension;
			Faulhaber(size_t _dimension);
			~Faulhaber()
			{
				delete [] polynoms;
			}
			const Polynom& operator()(size_t i) const;
			const static Faulhaber f;
		private:
			Polynom*const polynoms;

	};

}//namespace
#endif//guard

