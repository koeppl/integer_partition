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

#include "interval_partition.hpp"
#include <cassert>
#include <algorithm>
#include "binomial.hpp"
#include "faulhaber.hpp"
#include "sum_from_zero_to_upper.hpp"

using namespace IntervalPartition;

class Random
{
	gmp_randstate_t r_state;
	public:
	Random(unsigned long seed)
	{
		gmp_randinit_default(r_state);
		gmp_randseed_ui(r_state, seed);
	}
	~Random()
	{
		gmp_randclear(r_state);
	}
	Z get()
	{
		Z z;
		mpz_urandomb(z.get_mpz_t(),r_state,8);
		return z;
	}
};


Z fact(unsigned int n)
{
	 if(n < 2) return 1;
	 Z result(n);
	 while(n-- > 1) result *= n;
	 return result;
}
Q powerOf(Q value, unsigned int n)
{
	 if(n == 0) return 1;
	 Q result(value);
	 while(--n > 0) result *= value;
	 return result;
}
Q choose(unsigned int k, unsigned int j)
{
	 return Q(fact(k), fact(j) * fact(k-j));
}

/**
 * Evaluate the sum without using Faulhaber's formula,
 * e.g. iterating over the sum
 * Thus it's a barefoot approach without any technique
 */
Z barefoot(const IB& n, const unsigned int& p)
{
	Z res = 0;
	for(IB exp = 1; exp <= n; ++exp)
	{
		Z rop;
		mpz_pow_ui(rop.get_mpz_t(), exp.get_mpz_t(), p);
		res += rop;
	}
	return res;
}

Q sumFromZeroToUpperTest(const Polynom& p, const Z& upper)
{
	Q ret;
	for(Z k = 0; k <= upper; ++k)
		ret += p(k);
	return ret;
}

Q sumFromZMinusGammaToUpperTest(const Polynom& p, const Z& gamma, const Z& upper, const Z& z)
{
	Q ret;
	for(Z k = z-gamma; k <= upper; ++k)
		ret += p(k);
	return ret;
}



int main(void)
{
	Random random(1);
	
	{
		IntervalledPolynom ip;
		{
			Polynom p(5);
			p[0] = 1;
			p[1] = Q(3,2);
			p[2] = Q(1,2);
			ip.push_back(9, p);
		}
		{
			Polynom p(5);
			p[0] = -44;
			p[1] = 11;
			ip.push_back(30, p);
		}
		{
			Polynom p(5);
			p[0] = -479;
			p[1] = Q(81,2);
			p[2] = Q(-1,2);
			ip.push_back(40, p);
		}
		{
			Polynom p(5);
			p[0] = 341;
			ip.push_back(49, p);
		}
		{
			Polynom p(5);
			p[0] = -884;
			p[1] = Q(99,2);
			p[2] = Q(-1,2);
			ip.push_back(59, p);
		}
		{
			Polynom p(5);
			p[0] = 946;
			p[1] = -11;
			ip.push_back(80, p);
		}
		{
			Polynom p(5);
			p[0] = 4186;
			p[1] = Q(-183,2);
			p[2] = Q(1,2);
			ip.push_back(90, p);
		}
		{
			unsigned int dimensional_upper_bounds[] = {30, 50, 10};
			IntervalledPolynom intervalledPolynom = generateIntervalPartition(dimensional_upper_bounds, sizeof(dimensional_upper_bounds)/sizeof(unsigned int));
			assert(intervalledPolynom == ip);
		}
		{
			unsigned int dimensional_upper_bounds[] = {10, 30, 50};
			const size_t dimensional_upper_bounds_length = sizeof(dimensional_upper_bounds)/sizeof(unsigned int);
			const unsigned int maxdim = [&] () -> unsigned int
			{
				unsigned int maxdim_ = 0;
				for(const unsigned int& dim : dimensional_upper_bounds) maxdim_ += dim;
				return maxdim_;
			}();
			do
			{
				IntervalledPolynom intervalledPolynom = generateIntervalPartition(dimensional_upper_bounds, dimensional_upper_bounds_length);
				for(size_t x = 0; x <= maxdim; ++x)
					assert(intervalledPolynom(x) == ip(x));
			}
			while(std::next_permutation(dimensional_upper_bounds, dimensional_upper_bounds+dimensional_upper_bounds_length));
		}
	}

	// TESTS
	
	for(size_t x = 2; x < 100; ++x)
	for(size_t gamma = 0; gamma < x; ++gamma)
	for(size_t upper = x-gamma+1; upper < 100; ++upper)
	for(size_t p = 1; p < 20; ++p)
	{
		Polynom pol(p);
		for(size_t i = 0; i < p; ++i)
			pol[i] = random.get();
		Q res1 = sumFromZMinusGammaToUpper(pol,gamma,upper)(x);
		Q res2 = sumFromZMinusGammaToUpperTest(pol,gamma, upper, x);
//		std::cout << "SumToUpper: " << "x=" << x << ",gamma=" << gamma << "," << "upper=" << upper << ",p=" << p << "\t" << res1 << " " << res2 << std::endl;
		assert(res1 == res2);
	}	


	for(size_t i = 1; i < Binomial::b.dimension; ++i)
	for(size_t j = 0; j < i; ++j)
	{
		Q c = choose(i,j);
		Q b = Binomial::b(i,j);
		mpq_canonicalize(c.get_mpq_t());
		assert(mpq_equal(c.get_mpq_t(), b.get_mpq_t()));
	}	

	for(size_t x = 1; x < 100; ++x)
	for(size_t p = 0; p < 20; ++p)
	{
		Polynom pol = Faulhaber::f(p);
		Q res = pol(x);
		assert(res.get_den() == 1);
		Z resz = res.get_num();
		assert(resz == barefoot(x, p));
	}
	for(size_t x = 2; x < 100; ++x)
	for(size_t p = 1; p < 20; ++p)
	{
		Polynom pol(p);
		for(size_t i = 0; i < p; ++i)
			pol[i] = random.get();
		Q res1 = SumFromZeroToUpper::s(pol)(x);
		Q res2 = sumFromZeroToUpperTest(pol,x);
		assert(res1 == res2);
	}	

	for(size_t x = 2; x < 100; ++x)
	for(size_t gamma = 0; gamma < x; ++gamma)
	for(size_t p = 1; p < 20; ++p)
	{
		Polynom pol(p);
		for(size_t i = 0; i < p; ++i)
			pol[i] = random.get();
		Q res1 = sumFromZeroToZMinusGamma(pol,gamma)(x);
		Q res2 = sumFromZeroToUpperTest(pol,x-gamma);
		assert(res1 == res2);
	}	

	//END TESTS
	
	return 0;
}
