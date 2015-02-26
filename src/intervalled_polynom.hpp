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
 * @file intervalled_polynom.hpp
 * @brief A piecewise-defined polynomial function
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
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


	/** 
	 * A piecewise-defined polynomial function
	 * Each interval is represented by a single integer intervalbounds[i]
	 * The integers \f$ i_0, ..., i_{n-1} \f$ stored in intervalbounds span the intervals
	 * \f$ (-\infty, i_0], (i_0, i_1], ..., (i_{n-2}, i_{n-1}] $\f
	 */
	class IntervalledPolynom
	{
		private:
			vektor<IB> intervalbounds;
			vektor<Polynom> polynoms;
		public:

			/** 
			 * Returns the polynom that coincides with this polynomial at that given point.
			 * Takes \f$ O(\log n) \f$ time, where \f$ n \f$ is the number of different intervals.
			 * @pre The vector intervalbounds is sorted ascendingly.
			 * 
			 * @param point a point of the domain of this polynomial
			 * 
			 * @return Polynom defined at that point
			 */
			const Polynom& at(const IB& point) const;

			/** 
			 * Adds for a new interval a polynom on which this polynomial shall coincide with.
			 * The intervals are built from smallest to largest value.
			 * That means that you have to sort if you need to add intervalbounds in a non-asencending order!
			 * @pre intervalbound is an integer larger than any other integer that is used before any other call of this method.
			 * @post Both vectors intervalbounds and polynoms are appended by the new value (intervalbound, polynom)
			 * 
			 * @param intervalbound for an interval [a,b] the right value b
			 * @param polynom The polynom with which this polynomial shall coincide
			 */
			void push_back(const IB& intervalbound, const Polynom& polynom);
			/**
			 * Consumes additionally the polynom.
			 */
			void push_back(const IB& intervalbound, const Polynom&& polynom);

			/** 
			 * Swaps contents with another polynomial
			 */
			void swap(IntervalledPolynom& o);

			/** 
			 * Evaluates the polynomial at position x
			 * 
			 * @return The evaluated value
			 */
			Q operator()(const Z& x) const;
			bool operator==(IntervalledPolynom& o);

		friend std::ostream& operator<<(std::ostream& os, const IntervalledPolynom& ip);
	};

}//namespace
std::ostream& operator<<(std::ostream& os, const IntervalPartition::IntervalledPolynom& ip);
#endif//guard

