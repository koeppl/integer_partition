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
#include "sum_from_zero_to_upper.hpp"
#include "faulhaber.hpp"
#include "binomial.hpp"

namespace IntervalPartition
{

SumFromZeroToUpper SumFromZeroToUpper::s;
const Polynom& SumFromZeroToUpper::operator()(const Polynom& p)
{
	std::map<Polynom, Polynom>::const_iterator it = cache.find(p);
	if(it != cache.end()) return it->second;
	cache.insert(std::pair<Polynom,Polynom>(p, Polynom(p.size()+1)));
	Polynom& ret = cache.find(p)->second;
	
	for(size_t l = 0; l < p.size()+1; ++l) // l : Index of faulhaber's polynom
	{
		Q& koeff = ret[l];
		for(size_t j = 0; j < p.size(); ++j) // j : Index of p's coeffs
		{
			const Polynom& faulhaberpolynom = Faulhaber::f(j);
			if(faulhaberpolynom.size() <= l) continue;
			koeff += faulhaberpolynom.at(l) * p[j];
		}
	}
	ret[0] += p[0];
	return ret;
}

Polynom sumFromZeroToZMinusGamma(const Polynom& p, const Z& gamma)
{
	Polynom ret(p.size()+1);
	const Polynom& summedUp = SumFromZeroToUpper::s(p);
	for(size_t m = 0; m < p.size()+1; ++m) // m: index of leibniz binomial formula
	{
		Q& coeff = ret[m];
		Z gammapot = 1;
		for(size_t l = m; l < p.size()+1; ++l)
		{
			coeff += summedUp[l] * Binomial::b(l,m) * gammapot;
//			std::cout << "l=" << l << " m=" << m << " coeff=" << coeff << " gammapot=" << gammapot << std::endl;
	//		if( (l-m) % 2) coeff = -coeff;
			gammapot *= -gamma;
		}
	}
//	std::cout << "ret=" << ret << std::endl;
	return ret;
}

}//namespace
