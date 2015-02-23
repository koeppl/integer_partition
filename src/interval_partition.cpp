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
	 * The first term of the proof from Theorem 4.7
	 * It is computed by the complete sum to the upper bound minus the sum from $0$ to $z-\gamma$.
	 */
	Polynom sumFromZMinusGammaToUpper(const Polynom& p, const Z& gamma, const Z& upper)
	{
		Q highest = SumFromZeroToUpper::s(p)(upper);
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
		DEBUG_MSG("Interval Partitioning started");
#ifndef NDEBUG
		for(size_t i = 0; i < dimensions; ++i) {
			DCHECK_GT(dimensional_upper_bounds[i], 0);


		}
#endif
		/**
		 * Every dimension has to be > 0!
		*/

	//	vektor<std::pair<vektor<unsigned int>, vektor<unsigned int> > solution;
		vektor<IB> intervalbounds;

		intervalbounds.push_back(dimensional_upper_bounds[0]);

		IntervalledPolynom intervalledPolynom;
		{
			Polynom pol(1);
			pol[0] = 1;
			intervalledPolynom.push_back(dimensional_upper_bounds[0], pol);
		}//! This is exactly the induction base of Theorem 4.7

		for(size_t k = 1; k < dimensions; ++k)
		{
			DEBUG_MSG("k: " << k);
			assert(dimensional_upper_bounds[k] > 0);
			vektor<IB> help_intervalbounds;
			help_intervalbounds.push_back(dimensional_upper_bounds[k] - 1);

			for(const IB& old_intervalbound : intervalbounds)
				help_intervalbounds.push_back(old_intervalbound + dimensional_upper_bounds[k]);

			vektor<IB> tmp_intervalbounds; //! in this array the interval bounds of the next round (k+1) will be stored
			IntervalledPolynom tmp_intervalledPolynom; //! this will be the polynom of the next round (k+1)

			
			for(size_t i = 0, j = 0, witness_left = 0, witness_right = 0, last_upper_bound_index = 0, last_lower_bound_index = 0; j < help_intervalbounds.size();)
			{
				const IB intervalbound = i < intervalbounds.size() ?  intervalbounds[i] : 0;
				const IB help_intervalbound = help_intervalbounds[j];
				const IB last_upper_bound = last_upper_bound_index == 0 ? 0 : (last_upper_bound_index > intervalbounds.size() ? intervalbounds[last_upper_bound_index-2] : intervalbounds[last_upper_bound_index-1]);
				const IB last_lower_bound = last_lower_bound_index == 0 ? 0 : (last_lower_bound_index > intervalbounds.size() ? intervalbounds[last_lower_bound_index-2] : intervalbounds[last_lower_bound_index-1]);

				DEBUG_MSG("");
				DEBUG_MSG("k: " << k << ", i: " << i << ", j: " << j);
				if(i < intervalbounds.size())
				{
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
						DEBUG_MSG("Fall 1, addiere " << intervalbound);
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
						DEBUG_MSG("Fall 2, addiere " << help_intervalbound);
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
						DEBUG_MSG("Fall 3, addiere " << help_intervalbound);
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
					DEBUG_MSG("Fall 4, addiere " << help_intervalbound);
					// jetzt summe mit Anfangsintervall [...,
					// _old_intervals[j-1]]
					// und Endintervall [_old_intervals[i-1], ...]auswerten
					// Ergebnis zu _new_coefficients hinzufuegen
					++j;
				}
				DEBUG_MSG("intervalledPolynom: " << intervalledPolynom);
				DEBUG_MSG("tmp_intervalledPolynom: " << tmp_intervalledPolynom);
				DEBUG_MSG("Witness: " << "[" << witness_left << ", " << witness_right << "]");
				DEBUG_MSG("tmp_intervalbounds: " << tmp_intervalbounds);
				DEBUG_MSG("intervalbounds: " << intervalbounds);


				if(witness_left == witness_right)
				{
					if(witness_left <= intervalbounds.size() && witness_left > 0)
					{
						const IB& current_intervalbound = intervalbounds[witness_left-1]; //!
						DEBUG_MSG("Intervalbound: " << current_intervalbound);
						const Polynom& toSum = intervalledPolynom.at(current_intervalbound);
						const Polynom& upper = SumFromZeroToUpper::s(toSum);
						const Polynom lower = sumFromZeroToZMinusGamma(toSum, dimensional_upper_bounds[k]+1);
						Polynom together(std::max(upper.size(), lower.size())+1);
						for(size_t l = 0; l < together.size(); ++l)
						{
							if(l < upper.size()) together[l] += upper[l];
							if(l < lower.size()) together[l] -= lower[l];
						}
						DEBUG_MSG("Together Sum: " << together);
						tmp_intervalledPolynom.push_back(tmp_intervalbounds.back(), together);
					}
					else
						tmp_intervalledPolynom.push_back(tmp_intervalbounds.back(), Polynom::zero);

				}
				else
				{
					const Polynom lower_sum = (witness_left < 1 || witness_left-1 >= intervalbounds.size()) ? Polynom::zero : [&] () -> Polynom
					{
						const IB& lower_sum_intervalbound = intervalbounds[witness_left-1]; //!
						DEBUG_MSG("Lower Interval: " << lower_sum_intervalbound);
						return sumFromZMinusGammaToUpper(intervalledPolynom.at(lower_sum_intervalbound), dimensional_upper_bounds[k], lower_sum_intervalbound);
					}();
					DEBUG_MSG("Lower Sum: " << lower_sum);
					
					
					const Z const_sum = ( (witness_left > 0 || witness_right > 0) && witness_right-witness_left > 1) ? [&] () -> Z
					{
						Q const_sumq;
						for(size_t constinterval = witness_left; constinterval <= witness_right-2; ++constinterval)
						{
							if(constinterval >= intervalbounds.size()) break; //!
							const IB& intervalupper_bound = intervalbounds[constinterval]; //!
							const Polynom& summedUp = SumFromZeroToUpper::s(intervalledPolynom.at(intervalupper_bound));
							const_sumq += summedUp(intervalupper_bound);
							if(constinterval > 0)
							{
								const IB& intervallower_bound = intervalbounds[constinterval-1]; //!
								const_sumq -= summedUp(intervallower_bound);
							}
						}
						mpq_canonicalize(const_sumq.get_mpq_t());
						assert(const_sumq.get_den() == 1);
						return const_sumq.get_num();
					}() : 0;
					DEBUG_MSG("Constant Sum: " << const_sum);


					const Polynom upper_sum = (witness_right-1 >= intervalbounds.size()) ? Polynom::zero : [&] () -> Polynom
					{
						const IB& upper_sum_intervalupper_bound = intervalbounds[witness_right-1]; //!
						Polynom upper_sum_pol = SumFromZeroToUpper::s(intervalledPolynom.at(upper_sum_intervalupper_bound));
						if(witness_right > 1)
						{
							const IB& upper_sum_intervallower_bound = intervalbounds[witness_right-2];
							upper_sum_pol[0] -= upper_sum_pol(upper_sum_intervallower_bound);
						}
						DEBUG_MSG("Upper Interval: " << upper_sum_intervalupper_bound);
						return upper_sum_pol;
					}();
					DEBUG_MSG("Upper Sum: " << upper_sum);

					Polynom together(std::max(upper_sum.size(), lower_sum.size())+1);
					//TODO: unnecessary complicated!
					for(size_t l = 0; l < together.size(); ++l)
					{
						if(l < lower_sum.size()) together[l] += lower_sum[l];
						if(l < upper_sum.size()) together[l] += upper_sum[l];
					}
					together[0] += const_sum;
					DEBUG_MSG("Together Sum: " << together);
					tmp_intervalledPolynom.push_back(tmp_intervalbounds.back(), together);
				}

			}
			intervalbounds.swap(tmp_intervalbounds);
			intervalledPolynom.swap(tmp_intervalledPolynom);
			DEBUG_MSG("_old_intervals: " << intervalbounds);
		}
		DEBUG_MSG("Resulting Intervalled Polynom: " << intervalledPolynom);
		return intervalledPolynom;
	}

}//namespace
