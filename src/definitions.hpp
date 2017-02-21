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
 * @file definitions.hpp
 * @brief Definitions used across the library
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
 */

#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <vector>
#include <gmpxx.h>
#include <iostream>
#include "checked_vector.hpp"

#define BINOMIAL_DIM 200
#define BERNOULLI_DIM BINOMIAL_DIM-2
#define FAULHABER_DIM BINOMIAL_DIM-4

typedef mpq_class Q;
typedef mpz_class Z;
typedef mpz_class IB;
//template <typename T> using vektor = std::vector<T>;
template <typename T> using vektor = IntervalPartition::checked_vector<T>;

extern const Z Z_zero;

#ifndef NDEBUG
#include "debug.hpp"
#endif

#endif//guard
