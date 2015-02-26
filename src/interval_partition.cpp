/* Integer Partition
 * Computes the number of possible ordered integer partitions with upper bounds
 * Copyright (C) 2013 Dominik Köppl, Roland Glück
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
#include "sum_from_zero_to_upper.hpp"

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	os << "(";
	for(size_t i = 0; i < v.size(); ++i)
	{
		os << v[i] << " ";
	}
	os << ")";
	return os;
}

namespace IntervalPartition
{

	/** 
	 * Checks whether for a sequence (a_0,...,a_n) in container t_cont the condition
	 * \f$ p(a_{i+1}, a_i) \f$ holds for every \f$ 0 \le i \le n-1 \f$
	 *
	 * @param cont the container
	 * @param p the binary predicate, for example \class std::greater imposes a strictly increasing sequence.
	 *
	 * @tparam t_cont container class with forward iterator support
	 * @tparam BinaryPrecidate Predicate with signature \code bool pred(const t_cont::value_type &a, const t_cont::value_type &b); \endcode
	 * 
	 * @return true if the above condition holds
	 */
	template<class t_cont, class BinaryPredicate>
	bool has_ordering(const t_cont& cont, BinaryPredicate p) {
		typename t_cont::const_iterator it = cont.begin();
		typename t_cont::const_iterator it2 = it;
		while( ++it != cont.end()) {
			if(!p(*it, *(it2++) )) return false;
		}
		return true;
	}

	/** 
	 * The first term of the proof from Theorem 4.7
	 * It is computed by the complete sum to the upper bound minus the sum from $0$ to $z-\gamma$.
	 */
	Polynom sumFromZMinusGammaToUpper(const Polynom& p, const Z& gamma, const Z& upper) {
		const Q highest = SumFromZeroToUpper::s(p)(upper);
		Polynom diff = sumFromZeroToZMinusGamma(p,gamma+1);
		for(Q& coeff : diff)
			coeff = -coeff;
		diff[0] += highest;
		return diff;
	}

	/**
	 * Generates an intervalled polynom based on the interval bounds given as parameter
	 * @pre \code length(dimensional_upper_bounds) == dimensions \endcode has to hold.
	 */
	IntervalledPolynom generateIntervalPartition(const unsigned int* const dimensional_upper_bounds, const size_t dimensions)
	{
		LOG(INFO) << "Interval Partitioning started";
#ifndef NDEBUG
		for(size_t i = 0; i < dimensions; ++i) {
			DCHECK_GT(dimensional_upper_bounds[i], 0) << "Every dimensional upper bound has to be > 0";
		}
#endif
	//	vektor<std::pair<vektor<unsigned int>, vektor<unsigned int> > solution;
		vektor<IB> intervalbounds;

		intervalbounds.push_back(dimensional_upper_bounds[0]);

		IntervalledPolynom intervalledPolynom;
		{
			Polynom pol(1);
			pol[0] = 1;
			intervalledPolynom.push_back(dimensional_upper_bounds[0], pol);
		}//!< This is exactly the induction base of Theorem 4.7

		for(size_t k = 1; k < dimensions; ++k)
		{
			LOG(INFO) << "k: " << k;
			vektor<IB> help_intervalbounds;
			if(dimensional_upper_bounds[k] > 1) // we do now want a help_intervalbounds-value of 0
				help_intervalbounds.push_back(dimensional_upper_bounds[k] - 1);

			for(const IB& old_intervalbound : intervalbounds)
				help_intervalbounds.push_back(old_intervalbound + dimensional_upper_bounds[k]);

			vektor<IB> tmp_intervalbounds; //! in this array the interval bounds of the next round (k+1) will be stored
			IntervalledPolynom tmp_intervalledPolynom; //! this will be the polynom of the next round (k+1)

			
			/**
			 * i : intervalbounds[] index
			 * j : help_intervalbounds[] index
			 * last_upper_bound_index, last_lower_bound_index : intervalbounds[] index
			 * last_upper_bound, last_lower_bound : help_intervalbounds[] value
			 *
			 */
			for(size_t i = 0, j = 0, witness_left = 0, witness_right = 0, last_upper_bound_index = 0, last_lower_bound_index = 0; j < help_intervalbounds.size();)
			{
				//TODO This SUCKS:
				DCHECK_EQ(witness_right, last_upper_bound_index);
				DCHECK_EQ(witness_left, last_lower_bound_index);  

				const IB& intervalbound = i < intervalbounds.size() ?  intervalbounds[i] : 0;
				const IB& help_intervalbound = help_intervalbounds[j];
				//const IB& last_upper_bound = last_upper_bound_index == 0 ? 0 : (last_upper_bound_index > intervalbounds.size() ? intervalbounds[last_upper_bound_index-2] : intervalbounds[last_upper_bound_index-1]);
				//const IB& last_lower_bound = last_lower_bound_index == 0 ? 0 : (last_lower_bound_index > intervalbounds.size() ? intervalbounds[last_lower_bound_index-2] : intervalbounds[last_lower_bound_index-1]);
				const IB& last_upper_bound = last_upper_bound_index == 0 ? 0 : (last_upper_bound_index > intervalbounds.size() ? intervalbounds[intervalbounds.size()-1] : intervalbounds[last_upper_bound_index-1]);
				const IB& last_lower_bound = last_lower_bound_index == 0 ? 0 : (last_lower_bound_index > intervalbounds.size() ? intervalbounds[intervalbounds.size()-1] : intervalbounds[last_lower_bound_index-1]);
				DCHECK_LE(last_lower_bound_index, intervalbounds.size()+1);
				DCHECK_LE(last_upper_bound_index, intervalbounds.size()+1);

				LOG(INFO) << "k: " << k << ", i: " << i << ", j: " << j;
				if(i < intervalbounds.size())
				{
					/**
					 * We are examining intervalbounds[i] and help_intervalbounds[j] and pop
					 * that value which is smaller (popping: increment either i or j)
					 * If both values are the same, we increment both i and j.
					 * min(intervalbounds[i], help_intervalbounds[j]) is appended to tmp_intervalbounds
					 */
					if(intervalbound < help_intervalbound)
					{
						tmp_intervalbounds.push_back(intervalbound);
						if(last_upper_bound < intervalbound)
						{
							++witness_right;
							++last_upper_bound_index;
						}
						if(intervalbound > last_lower_bound + dimensional_upper_bounds[k] ) //TODO: with else?
						{
							++witness_left;
							++last_lower_bound_index;
						}
						LOG(INFO) << "Fall 1, addiere " << intervalbound;
						// jetzt summe mit Anfangsintervall [...,
						// _old_intervals[j-1]]
						// und Endintervall [_old_intervals[i-1], ...]auswerten
						// Ergebnis zu _new_coefficients hinzufuegen
						++i;
					}
					else if(intervalbound > help_intervalbound)
					{
						tmp_intervalbounds.push_back(help_intervalbound);
						if(last_upper_bound < help_intervalbound)
						{
							++witness_right;
							++last_upper_bound_index;
						}
						if(help_intervalbound > last_lower_bound + dimensional_upper_bounds[k])
						{
							++witness_left;
							++last_lower_bound_index;
						}
						LOG(INFO) << "Fall 2, addiere " << help_intervalbound;
						// jetzt summe mit Anfangsintervall [...,
						// _old_intervals[j-1]]
						// und Endintervall [_old_intervals[i-1],
						// ...]auswerten
						// Ergebnis zu _new_coefficients hinzufuegen
						++j;
					}
					else
					{
						tmp_intervalbounds.push_back(intervalbound);
						if(intervalbound > last_upper_bound)
						{
							++witness_right;
							++last_upper_bound_index;
						}
						if(intervalbound > last_lower_bound + dimensional_upper_bounds[k])
						{
							++witness_left;
							++last_lower_bound_index;
						}
						LOG(INFO) << "Fall 3, addiere " << help_intervalbound;
						// jetzt summe mit Anfangsintervall [...,
						// _old_intervals[j-1]]
						// und Endintervall [_old_intervals[i-1],
						// ...]auswerten
						// Ergebnis zu _new_coefficients hinzufuegen
						++i;
						++j;
					}
				} 
				else 
					/** We have already taken every element of intervalbounds[]. Because help_intervalbounds[] has some larger values, these have to 
					 *  be examined:
					 */
				{
					tmp_intervalbounds.push_back(help_intervalbound);
					if(help_intervalbound > last_upper_bound && last_upper_bound_index <= intervalbounds.size())
					{
						++witness_right;
						++last_upper_bound_index;
					}
					if(help_intervalbound > last_lower_bound + dimensional_upper_bounds[k])
					{
						++witness_left;
						++last_lower_bound_index;
					}
					/*

					if(help_intervalbound <= static_cast<long>(intervalbounds[intervalbounds.size()-1]))
					{
						witness_right = intervalbounds.size();
						++witness_left;
					}
					else
					{
						if(witness_right == intervalbounds.size())
							++witness_right;
						else
							++witness_left;
					}
					*/
					LOG(INFO) << "Fall 4, addiere " << help_intervalbound;
					// jetzt summe mit Anfangsintervall [...,
					// _old_intervals[j-1]]
					// und Endintervall [_old_intervals[i-1], ...]auswerten
					// Ergebnis zu _new_coefficients hinzufuegen
					++j;
				}
				LOG(INFO) << "intervalledPolynom: " << intervalledPolynom;
				LOG(INFO) << "tmp_intervalledPolynom: " << tmp_intervalledPolynom;
				LOG(INFO) << "Witness: " << "[" << witness_left << ", " << witness_right << "]";
				LOG(INFO) << "tmp_intervalbounds: " << tmp_intervalbounds;
				LOG(INFO) << "intervalbounds: " << intervalbounds;
				
				// TODO: Make the lower part a function for itself

				/**
				 * Here, we compute the new polynomial over the bounds 
				 * intervalbounds[witness_left-1] to intervalbounds[witness_right-1]
				 */

				if(witness_left == witness_right)
				{
					/**
					 * If witness_left == witness_right, then there is no sum over a constant range.
					 * We calculate
					 * \f$ \sum_{k= z - i_n}^z p(k) = \sum_{k=0}^z p(k) - \sum_{k=0}^{z-i_n-1} p(k) \f$
					 * with i_n = witness_left-2 and p = intevalledPolynom.at(intervalbounds[i_n+1])
					 *
					 */
					if(witness_left <= intervalbounds.size() && witness_left > 0)
					{
						const IB& current_intervalbound = intervalbounds[witness_left-1]; //!
						LOG(INFO) << "Intervalbound: " << current_intervalbound;
						const Polynom& toSum = intervalledPolynom.at(current_intervalbound);
						const Polynom& upper = SumFromZeroToUpper::s(toSum);
						const Polynom lower = sumFromZeroToZMinusGamma(toSum, dimensional_upper_bounds[k]+1);
						{// adds upper + lower to tmp_intervalledPolynom // TODO: make operator+ out of it
							Polynom together(std::max(upper.size(), lower.size())+1);
							for(size_t l = 0; l < together.size(); ++l)
							{
								if(l < upper.size()) together[l] += upper[l];
								if(l < lower.size()) together[l] -= lower[l];
							}
							LOG(INFO) << "Together Sum: " << together;
							tmp_intervalledPolynom.push_back(tmp_intervalbounds.back(), std::move(together));
						}
					}
					else 
						tmp_intervalledPolynom.push_back(tmp_intervalbounds.back(), Polynom::zero);

					DCHECK_GE(tmp_intervalledPolynom.at(tmp_intervalbounds.back())(tmp_intervalbounds.back()), 0); // Invariant: polynomial is non-negative
				}
				else
				{
					/**
					 * Left-Witness-Intervall: intervalbounds[witness_left-1]
					 * Left-Witness: [s,t]
					 *
					 * \f$ \sum_{l= z - i_k }^{t} p(k) \f$
					 * where p is intervallPolynom.at(Left-Witness)
					 */
					const Polynom lower_sum = (witness_left < 1 || witness_left-1 >= intervalbounds.size()) ? Polynom::zero : [&] () -> Polynom
					{
						const IB& lower_sum_intervalbound = intervalbounds[witness_left-1]; //!
						LOG(INFO) << "Lower Interval: " << lower_sum_intervalbound;
						return sumFromZMinusGammaToUpper(intervalledPolynom.at(lower_sum_intervalbound), dimensional_upper_bounds[k], lower_sum_intervalbound);
					}();
					LOG(INFO) << "Lower Sum: " << lower_sum;
					DCHECK_GE(lower_sum(tmp_intervalbounds.back()), 0) << "Lower Sum is Negative!"; // Invariant: polynomial is non-negative
					
					/**
					 * Let:
					 * l: witness_left  -1 (lower bound)
					 * u: witness_right -1 (upper bound)
					 *
					 * Calculate constant sum of intervalledPolynom.at(k) over \f$ k \in I^n_j \f$ with \f$ l < j < u \f$
					 * This range is split up by \f$ \sum_k={s_j}^{t_j} p(x) \f$ for every \f$ l < j < u\$ and $I^n_j = (s_j, t_j) \f$
					 */
					const Z const_sum = ( (witness_left > 0 || witness_right > 0) && witness_right-witness_left > 1) ? [&] () -> Z
					{
						Q const_sumq;
						for(size_t constinterval = witness_left; constinterval <= witness_right-2; ++constinterval)
						{
							if(constinterval >= intervalbounds.size()) break; //!
							const IB& intervalupper_bound = intervalbounds[constinterval]; //!
							const Polynom& summedUp = SumFromZeroToUpper::s(intervalledPolynom.at(intervalupper_bound));
							const_sumq += summedUp(intervalupper_bound);
							if(constinterval > 0) //!< if we examine \f$I^n_j\f$ for a j > 0, then we have summed up to much -> subtract interval [0, s_j[
							{
								const IB& intervallower_bound = intervalbounds[constinterval-1]; //!
								const_sumq -= summedUp(intervallower_bound);
							}
						}
						mpq_canonicalize(const_sumq.get_mpq_t());
						DCHECK_EQ(const_sumq.get_den(), 1);
						return const_sumq.get_num();
					}() : 0;
					LOG(INFO) << "Constant Sum: " << const_sum;
					DCHECK_GE(const_sum, 0) << "Constant Sum is Negative!"; // Invariant: polynomial is non-negative


					/**
					 * Right-Witness-Intervall: intervalbounds[witness_right-1]
					 * Right-Witness: [s,t]
					 *
					 * \f$ \sum_{k= s}^{z} p(k) \f$
					 * where p is intervallPolynom.at(Right-Witness)
					 */
					const Polynom upper_sum = (witness_right-1 >= intervalbounds.size()) ? Polynom::zero : [&] () -> Polynom
					{
						const IB& upper_sum_intervalupper_bound = intervalbounds[witness_right-1]; //!
						LOG(INFO) << "Upper Interval: " << upper_sum_intervalupper_bound;
						Polynom upper_sum_pol = SumFromZeroToUpper::s(intervalledPolynom.at(upper_sum_intervalupper_bound));
						LOG(INFO) << "upper_sum_pol pre " << upper_sum_pol;
						if(witness_right > 1)
						{
							const IB& upper_sum_intervallower_bound = intervalbounds[witness_right-2];
							upper_sum_pol[0] -= upper_sum_pol(upper_sum_intervallower_bound);
						}
						LOG(INFO) << "upper_sum_pol post " << upper_sum_pol;
						return upper_sum_pol;
					}();
					LOG(INFO) << "Upper Sum: " << upper_sum;
					DCHECK_GE(upper_sum(tmp_intervalbounds.back()), 0) << "Upper Sum is Negative!"; // Invariant: polynomial is non-negative

					Polynom together(std::max(upper_sum.size(), lower_sum.size())+1);
					//TODO: unnecessary complicated!
					for(size_t l = 0; l < together.size(); ++l)
					{
						if(l < lower_sum.size()) together[l] += lower_sum[l];
						if(l < upper_sum.size()) together[l] += upper_sum[l];
					}
					together[0] += const_sum;
					LOG(INFO) << "Together Sum: " << together;
					tmp_intervalledPolynom.push_back(tmp_intervalbounds.back(), std::move(together));
					DCHECK_GE(tmp_intervalledPolynom.at(tmp_intervalbounds.back())(tmp_intervalbounds.back()), 0); // Invariant: polynomial is non-negative
				}

			}

#ifndef NDEBUG
			DCHECK(has_ordering(tmp_intervalbounds, std::greater<IB>())); // Invariant: the numbers of tmp_intervalbounds are strict ascendending
			for(const auto& ibound : tmp_intervalbounds) { //Invariant: the piecewise-defined polynomial is non-negative.
				DCHECK_GE(tmp_intervalledPolynom.at(ibound)(ibound), 0);
			}
#endif

			intervalbounds.swap(tmp_intervalbounds);
			intervalledPolynom.swap(tmp_intervalledPolynom);
			LOG(INFO) << "_old_intervals: " << intervalbounds;
		}
		LOG(INFO) << "Resulting Intervalled Polynom: " << intervalledPolynom;
		return intervalledPolynom;
	}



}//namespace
