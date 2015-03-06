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
 * @file polynom.hpp
 * @brief Polynom Functional
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
 */
#ifndef POLYNOM_HPP
#define POLYNOM_HPP
#include "definitions.hpp"

namespace IntervalPartition
{

	/** 
	 * Uses std::vector for storing the coefficients.
	 * Overloads operator() for evaluationg the polynom at any point.
	 */
	class Polynom : public vektor<Q>
	{
		public:

			/** 
			 * Creates an empty polynom, i.e., all coefficienst are set to zero.
			 * @param psize Creates a polynom with psize coefficients
			 */
		Polynom(size_t psize) : vektor<Q>(psize) { }
		Polynom(size_t psize, const size_t fill) : vektor<Q>(psize) { Q f = fill; for(auto& i : *this) i = f;  }
		explicit Polynom(const Polynom& pol) : vektor<Q>(pol) {}
		Polynom(Polynom&& pol) : vektor<Q>(pol) {}
		Polynom() : vektor<Q>() {}

		/**
		 * Evaluates the polynom at position x
		 * Uses Horner method
		 * a_0 + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1}
		 * where std::vector stores the sequence (a_0, ..., a_{n-1})
		 * 
		 */
		Q operator()(const Q& x) const
		{
			Q res = at(size()-1);
			for(size_t i = 1; i < size(); ++i)
				res = res*x + at(size()-i-1);
			mpq_canonicalize(res.get_mpq_t());
			return res;
		}
		/**
		 * Canonicalizes the coefficients of the polynom.
		 * These could be unnecessarily great, because gnump does not 
		 * refactor rational numbers automatically.
		 */
		void canonicalize()
		{
			while(back() == 0) pop_back();
			for(Q& coeff : *this) mpq_canonicalize(coeff.get_mpq_t());
		}
		const static Polynom zero; //!< the polynom with 0 as coefficient
		const static Polynom one; //!< the polynom with 0 as coefficient
	};
	std::ostream& operator<<(std::ostream& os, const Polynom& v);

	Polynom operator-(const Polynom& a, const Polynom& b);
	Polynom operator+(const Polynom& a, const Polynom& b);
}//ns
#endif//guard

