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
#ifndef INTERVALL_PARTITION
#define INTERVALL_PARTITION
#include "intervalled_polynom.hpp"

namespace IntervalPartition
{
	IntervalledPolynom generateIntervalPartition(const unsigned int* const dimensional_upper_bounds, const size_t dimensions);
	Polynom sumFromZMinusGammaToUpper(const Polynom& p, const Z& gamma, const Z& upper);
}
#endif//guard

