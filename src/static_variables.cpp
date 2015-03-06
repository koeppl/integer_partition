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
#include "bernoulli.hpp"
#include "faulhaber.hpp"

/** Do NOT change the order of the variable intitialization below! **/
namespace IntervalPartition {


	const Binomial Binomial::b(BINOMIAL_DIM);

	static_assert(BERNOULLI_DIM+1 < BINOMIAL_DIM, "Bernoulli dimension must be less than Binomial dimension -1");
	const Bernoulli Bernoulli::b(Binomial::b, BERNOULLI_DIM);

	static_assert(FAULHABER_DIM+3 < BINOMIAL_DIM, "Faulhaber dimension must be less than Binomial dimension -3");
	const Faulhaber Faulhaber::f(FAULHABER_DIM);

}
