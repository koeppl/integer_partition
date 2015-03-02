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

#include "definitions.hpp"
#include "random.hpp"
#include <algorithm>
#include "binomial.hpp"
#include "faulhaber.hpp"
#include "sum_from_zero_to_upper.hpp"
#include <gtest/gtest.h>

using namespace IntervalPartition;


/**
 * Calculates n!
 */
Z fact(unsigned int n)
{
	 if(n < 2) return 1;
	 Z result(n);
	 while(n-- > 1) result *= n;
	 return result;
}
/**
 * Calculates value^n
 */
Q powerOf(Q value, unsigned int n)
{
	 if(n == 0) return 1;
	 Q result(value);
	 while(--n > 0) result *= value;
	 return result;
}
/**
 * Calculates \f$k \choose j\f$
 */
Q choose(unsigned int k, unsigned int j)
{
	 return Q(fact(k), fact(j) * fact(k-j));
}

/**
 * Evaluate the sum without using Faulhaber's formula,
 * e.g. iterating over the sum
 * Thus it's a naive approach without any technique
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
/**
 * Calculates \f$\sum_{k=0}^{\text{upper}} p(k)\f$
 */
Q sumFromZeroToUpperTest(const Polynom& p, const Z& upper)
{
	Q ret;
	for(Z k = 0; k <= upper; ++k)
		ret += p(k);
	return ret;
}

/**
 * Calculates \f$\sum_{k=z-\gamma}^{\text{upper}} p(k)\f$ naively
 *
 */
Q sumFromZMinusGammaToUpperTest(const Polynom& p, const Z& gamma, const Z& upper, const Z& z)
{
	Q ret;
	for(Z k = z-gamma; k <= upper; ++k)
		ret += p(k);
	return ret;
}



unsigned long seed = 0; //TODO




// TESTS
	

TEST(Polynom, SumFormZMinusGammaToUpper) {
	IntervalPartition::GMPRandom random(seed);
	for(size_t x = 2; x < 10; ++x)
	for(size_t gamma = 0; gamma < x; ++gamma)
	for(size_t upper = x-gamma+1; upper < 100; ++upper)
	for(size_t p = 1; p < 20; ++p)
	{
		Polynom pol(p);
		for(size_t i = 0; i < p; ++i)
			pol[i] = random.get();
		const Q res1 = sumFromZMinusGammaToUpper(pol,gamma,upper)(x);
		const Q res2 = sumFromZMinusGammaToUpperTest(pol,gamma, upper, x);
		ASSERT_EQ(res1, res2) << "SumToUpper: " << "x=" << x << ",gamma=" 
			<< gamma << "," << "upper=" << upper << ",p=" << p << "\t" << res1 << " " << res2 << std::endl;
	}	
}

TEST(Binomial, NaiveTest) {
	/** Testing \class Binomial
	 */
	for(size_t i = 1; i < Binomial::b.dimension; ++i)
	for(size_t j = 0; j < i; ++j) {
		Q c = choose(i,j);
		Q b = Binomial::b(i,j);
		mpq_canonicalize(c.get_mpq_t());
		ASSERT_TRUE(mpq_equal(c.get_mpq_t(), b.get_mpq_t()));
	}	
}
TEST(Polynom, Faulhaber) {
	for(size_t x = 1; x < 100; ++x)
	for(size_t p = 0; p < 20; ++p) {
		const Polynom& pol = Faulhaber::f(p);
		Q res = pol(x);
		ASSERT_EQ(res.get_den(), 1);
		Z resz = res.get_num();
		ASSERT_EQ(resz, barefoot(x, p));
	}
}
TEST(Polynom, SumFromZeroToUpper) {
	IntervalPartition::GMPRandom random(seed);
	for(size_t x = 2; x < 100; ++x)
	for(size_t p = 1; p < 20; ++p) {
		Polynom pol(p);
		for(size_t i = 0; i < p; ++i)
			pol[i] = random.get();
		const Q res1 = SumFromZeroToUpper::s(pol)(x);
		const Q res2 = sumFromZeroToUpperTest(pol,x);
		ASSERT_EQ(res1, res2);
	}	
}

TEST(Polynom, sumFromZeroToZMinusGamma) {
	IntervalPartition::GMPRandom random(seed);
	for(size_t x = 2; x < 100; ++x)
	for(size_t gamma = 0; gamma < x; ++gamma)
	for(size_t p = 1; p < 20; ++p) {
		Polynom pol(p);
		for(size_t i = 0; i < p; ++i)
			pol[i] = random.get();
		const Q res1 = sumFromZeroToZMinusGamma(pol,gamma)(x);
		const Q res2 = sumFromZeroToUpperTest(pol,x-gamma);
		ASSERT_EQ(res1, res2);
	}	
}

