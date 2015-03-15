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
#include "util.hpp"
#include "intervalled_polynom.hpp"
#include "sum_from_zero_to_upper.hpp"
#include "binomial.hpp"
#include <future>
#include <iterator>
#include <queue>
#include "sum_from_zero_threads.hpp"

namespace IntervalPartition {

	class PiecedPolyAsync {
		public:
		std::mutex evalMutex;
		PiecedPolyAsync(const PiecedPolyAsync&) = delete; // This class is to experimantal to be used for copying...
		PiecedPolyAsync() = default;
		Polynom& evaluate(size_t distance) {

			if(evaluated_polynoms[distance] == nullptr) {
				std::lock_guard<std::mutex> guard(evalMutex);
				if(evaluated_polynoms[distance] != nullptr) return *evaluated_polynoms[distance];
				try {
					evaluated_polynoms[distance] = polynoms[distance]->get();
					DCHECK(evaluated_polynoms[distance] != nullptr);// ((void*)evaluated_polynoms[distance]) << " is null!";
				} catch (const std::future_error& e) {
					DCHECK_EQ(false,true);
				}

			}
			return *evaluated_polynoms[distance];
		}
			~PiecedPolyAsync() {
				for(auto &p : evaluated_polynoms) if(p != nullptr) delete p;
				for(auto &p : polynoms) if(p != nullptr) delete p;
			}
			vektor<IB> intervalbounds;
			vektor<std::future<Polynom*>*> polynoms;
			vektor<Polynom*> evaluated_polynoms;
			const vektor<IB>& bounds() const{ return intervalbounds; }
			const Polynom& at(const IB& point);
			void push_back(const IB& intervalbound, std::future<Polynom*>&& polynom);
			void push_back(const IB& intervalbound, Polynom* polynom);
			void swap(PiecedPolyAsync& o);
			Q operator()(const Z& x);
			IntervalledPolynom extract() {
				IntervalledPolynom o;
				for(size_t i = 0; i < evaluated_polynoms.size(); ++i) {
					o.push_back( intervalbounds[i], std::move(evaluate(i)));
				}
				return o;
			}
	};

std::ostream& operator<<(std::ostream& os, PiecedPolyAsync& ip) {
	os << "IP: ";
	for(size_t i = 0; i < ip.intervalbounds.size(); ++i) {
		os << "{" << ip.intervalbounds[i] << " => " << ip.at(i) << "}";
		if(i+1 < ip.intervalbounds.size()) os << ", ";
	}
	return os;
}

const Polynom& PiecedPolyAsync::at(const IB& point) {
	if(point < 0) return Polynom::zero;
	LOG(INFO) <<  "Search " << point << " in " << intervalbounds;
	vektor<IB>::const_iterator it = lower_bound(intervalbounds.begin(), intervalbounds.end(), point);
	if(it != intervalbounds.end() && *it >= point) {
		const size_t distance = std::distance<vektor<IB>::const_iterator>(intervalbounds.begin(),it);
		Polynom& p = evaluate(distance);
		LOG(INFO) << "Found " << intervalbounds[distance] << "->" << p;
		return p;
	}
	LOG(INFO) << "Search failed !";
	return Polynom::zero;
}

Q PiecedPolyAsync::operator()(const Z& x) {
	return at(x)(x);
}
void PiecedPolyAsync::push_back(const IB& intervalbound, std::future<Polynom*>&& polynom) {
	intervalbounds.push_back(intervalbound);
	polynoms.push_back(new std::future<Polynom*>(std::move(polynom)));
	evaluated_polynoms.push_back(nullptr);
}
void PiecedPolyAsync::push_back(const IB& intervalbound, Polynom* polynom) {
	intervalbounds.push_back(intervalbound);
	polynoms.push_back(nullptr);
	evaluated_polynoms.push_back(polynom);
}
void PiecedPolyAsync::swap(PiecedPolyAsync& o) {
	intervalbounds.swap(o.intervalbounds);
	polynoms.swap(o.polynoms);
	evaluated_polynoms.swap(o.evaluated_polynoms);
}

	Polynom* sumParallelPolynomialOverWitnesses
	( const unsigned int& dimensional_upper_bound // dimensional_upper_bounds[k]
	, const vektor<IB>& intervalbounds 
	, const std::function<const Polynom&(const Polynom&)>& SumFromZeroToUpper
	, PiecedPolyAsync& intervalledPolynom
	, const size_t& witness_left_index
	, const size_t& witness_right_index
	)
		{
				/**
				 * Here, we compute the new polynomial over the bounds 
				 * intervalbounds[witness_left_index-1] to intervalbounds[witness_right_index-1]
				 */

				if(witness_left_index == witness_right_index) {
					/**
					 * If witness_left_index == witness_right_index, then there is no sum over a constant range.
					 * We calculate
					 * \f$ \sum_{k= z - i_n}^z p(k) = \sum_{k=0}^z p(k) - \sum_{k=0}^{z-i_n-1} p(k) \f$
					 * with i_n = witness_left_index-2 and p = intevalledPolynom.at(intervalbounds[i_n+1])
					 *
					 */
					DCHECK_LE(witness_left_index, intervalbounds.size());
					DCHECK_GT(witness_left_index, 0);
					const IB& current_intervalbound = intervalbounds[witness_left_index-1]; //!
					LOG(INFO) << "Intervalbound: " << current_intervalbound;
					const Polynom& toSum = intervalledPolynom.at(current_intervalbound);
					const Polynom& upper = SumFromZeroToUpper(toSum);
					const Polynom lower = sumFromZeroToZMinusGamma(toSum, dimensional_upper_bound+1, SumFromZeroToUpper, Binomial::b);
					Polynom* together = new Polynom(upper - lower);
					LOG(INFO) << "Together Sum: " << (*together);
					return together;
				} else {
					/**
					 * Left-Witness-Intervall: intervalbounds[witness_left_index-1]
					 * Left-Witness: [s,t]
					 *
					 * \f$ \sum_{l= z - i_k }^{t} p(k) \f$
					 * where p is intervallPolynom.at(Left-Witness)
					 */
					const Polynom lower_sum([&] () -> Polynom
					{
						if(witness_left_index < 1 || witness_left_index-1 >= intervalbounds.size()) {
							Polynom lower_sum_pol(Polynom::zero);
							return lower_sum_pol;
						}
						const IB& lower_sum_intervalbound = intervalbounds[witness_left_index-1]; //!
						LOG(INFO) << "Lower Interval: " << lower_sum_intervalbound;
						Polynom lower_sum_pol = sumFromZMinusGammaToUpper(intervalledPolynom.at(lower_sum_intervalbound), dimensional_upper_bound, lower_sum_intervalbound, SumFromZeroToUpper, Binomial::b);
						return lower_sum_pol;
					}());
					LOG(INFO) << "Lower Sum: " << lower_sum;
					
					/**
					 * Let:
					 * l: witness_left_index  -1 (lower bound)
					 * u: witness_right_index -1 (upper bound)
					 *
					 * Calculate constant sum of intervalledPolynom.at(k) over \f$ k \in I^n_j \f$ with \f$ l < j < u \f$
					 * This range is split up by \f$ \sum_k={s_j}^{t_j} p(x) \f$ for every \f$ l < j < u\$ and $I^n_j = (s_j, t_j) \f$
					 */
					const Z const_sum = ( (witness_left_index > 0 || witness_right_index > 0) && witness_right_index-witness_left_index > 1) ? [&] () -> Z
					{
						Q const_sumq;
						for(size_t constinterval = witness_left_index; constinterval <= witness_right_index-2; ++constinterval)
						{
							if(constinterval >= intervalbounds.size()) break; //!
							const IB& intervalupper_bound = intervalbounds[constinterval]; //!
							const Polynom& summedUp = SumFromZeroToUpper(intervalledPolynom.at(intervalupper_bound));
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
					 * Right-Witness-Intervall: intervalbounds[witness_right_index-1]
					 * Right-Witness: [s,t]
					 *
					 * \f$ \sum_{k= s}^{z} p(k) \f$
					 * where p is intervallPolynom.at(Right-Witness)
					 */
					const Polynom upper_sum( [&] () -> Polynom
					{
						if(witness_right_index-1 >= intervalbounds.size()) {
							Polynom upper_sum_pol(Polynom::zero);
							return upper_sum_pol;
						}
						const IB& upper_sum_intervalupper_bound = intervalbounds[witness_right_index-1]; //!
						LOG(INFO) << "Upper Interval: " << upper_sum_intervalupper_bound;
						Polynom upper_sum_pol = Polynom(SumFromZeroToUpper(intervalledPolynom.at(upper_sum_intervalupper_bound)));
						LOG(INFO) << "upper_sum_pol pre " << upper_sum_pol;
						if(witness_right_index > 1)
						{
							const IB& upper_sum_intervallower_bound = intervalbounds[witness_right_index-2];
							upper_sum_pol[0] -= upper_sum_pol(upper_sum_intervallower_bound);
						}
						LOG(INFO) << "upper_sum_pol post " << upper_sum_pol;
						DCHECK_GE(upper_sum_pol(intervalbounds[witness_right_index-1]), 0) << "Upper Sum is Negative!"; // Invariant: polynomial is non-negative
						return upper_sum_pol;
					}());
					LOG(INFO) << "Upper Sum: " << upper_sum;

					
					Polynom* together = new Polynom(upper_sum + lower_sum);
					(*together)[0] += const_sum;
					LOG(INFO) << "Together Sum: " << (*together);
					return together;
				}
		}



		struct Job {
			const size_t m_dimension;
			const size_t m_witness_left_index;
			const size_t m_witness_right_index;
			std::promise<Polynom*> m_result;
			Job(const size_t dimension, const size_t witness_left_index, const size_t witness_right_index, std::promise<Polynom*>&& result) 
				: m_dimension(dimension), m_witness_left_index(witness_left_index), m_witness_right_index(witness_right_index), m_result(std::move(result))
			{
			}
		};


	/**
	 * Generates an intervalled polynom based on the interval bounds given as parameter
	 * @pre \code length(dimensional_upper_bounds) == dimensions \endcode has to hold.
	 */
	IntervalledPolynom generateParallelIntervalPartition(const unsigned int* const dimensional_upper_bounds, const size_t dimensions, bool useSymmetry)
	{
		LOG(INFO) << "Interval Partitioning started";
#ifndef NDEBUG
		for(size_t i = 0; i < dimensions; ++i) {
			DCHECK_GT(dimensional_upper_bounds[i], 0) << "Every dimensional upper bound has to be > 0";
		}
#endif
		// maxdim is only used if useSymmetry to decide the cut-off
		const size_t maxdim = std::accumulate(dimensional_upper_bounds, dimensional_upper_bounds+dimensions,static_cast<size_t>(0)); 
		vektor<IB> intervalbounds;
		SumFromZeroCacherThreadSafe sumcacher; // instead of calling sumFromZeroToUpper each time, we cache its sum.
		std::function<const Polynom&(const Polynom&)> SumFromZeroToUpper = [&sumcacher] (const Polynom& a) -> const Polynom& { return sumcacher(a);};

		intervalbounds.push_back(dimensional_upper_bounds[0]);

		std::queue<Job> jobs;

		vektor<PiecedPolyAsync> piecewisePolynoms(dimensions);

		{
			Polynom* pol = new Polynom(1);
			(*pol)[0] = 1;
			piecewisePolynoms[0].push_back(dimensional_upper_bounds[0], pol);
		}//!< This is exactly the induction base of Theorem 4.7

		for(size_t k = 1; k < dimensions; ++k)
		{
			LOG(INFO) << "k: " << k;

			vektor<IB> tmp_intervalbounds; //! in this array the interval bounds of the next round (k+1) will be stored
//			PiecedPolyAsync tmp_intervalledPolynom; //! this will be the polynom of the next round (k+1)
			vektor<IB> help_intervalbounds;
			const unsigned int& dimensional_upper_bound = dimensional_upper_bounds[k];

			//if(dimensional_upper_bound > 1) 
			{ // we do now want a help_intervalbounds-value of 0
				help_intervalbounds.push_back(dimensional_upper_bound - 1);
			}


			for(const IB& old_intervalbound : intervalbounds)
				help_intervalbounds.push_back(old_intervalbound + dimensional_upper_bound);


			
			/**
			 * i : intervalbounds[] index
			 * j : help_intervalbounds[] index
			 * witness_right_index, witness_left_index : intervalbounds[] index
			 * witness_right, witness_left : help_intervalbounds[] value
			 *
			 */
			for(size_t i = 0, j = 0, witness_left_index = 0, witness_right_index = 0; j < help_intervalbounds.size();)
			{

				const IB& intervalbound = i < intervalbounds.size() ?  intervalbounds[i] : Z_zero;
				const IB& help_intervalbound = help_intervalbounds[j];
				
				const IB& witness_right = get_witness(witness_right_index, intervalbounds);
				const IB& witness_left = get_witness(witness_left_index, intervalbounds);
				if(useSymmetry && witness_left >= maxdim/2) break;
				DCHECK_LE(witness_left_index, intervalbounds.size()+1);
				DCHECK_LE(witness_right_index, intervalbounds.size()+1);

				LOG(INFO) << "k: " << k << ", i: " << i << ", j: " << j;
				if(i < intervalbounds.size()) {
					/**
					 * We are examining intervalbounds[i] and help_intervalbounds[j] and pop
					 * that value which is smaller (popping: increment either i or j)
					 * If both values are the same, we increment both i and j.
					 * min(intervalbounds[i], help_intervalbounds[j]) is appended to tmp_intervalbounds
					 */
					if(intervalbound < help_intervalbound) {
						tmp_intervalbounds.push_back(intervalbound);
						if(witness_right < intervalbound) {
							++witness_right_index;
						}
						if(intervalbound > witness_left + dimensional_upper_bound ) {
							++witness_left_index;
						}
						LOG(INFO) << "Fall 1, addiere " << intervalbound;
						++i;
					}
					else if(intervalbound > help_intervalbound)
					{
						tmp_intervalbounds.push_back(help_intervalbound);
						if(witness_right < help_intervalbound) {
							++witness_right_index;
						}
						if(help_intervalbound > witness_left + dimensional_upper_bound) {
							++witness_left_index;
						}
						LOG(INFO) << "Fall 2, addiere " << help_intervalbound;
						++j;
					}
					else
					{
						tmp_intervalbounds.push_back(intervalbound);
						if(intervalbound > witness_right) {
							++witness_right_index;
						}
						if(intervalbound > witness_left + dimensional_upper_bound) {
							++witness_left_index;
						}
						LOG(INFO) << "Fall 3, addiere " << intervalbound;
						++i;
						++j;
					}
				} else {
					/** We have already taken every element of intervalbounds[]. Because help_intervalbounds[] has some larger values, these have to 
					 *  be examined:
					 */
					tmp_intervalbounds.push_back(help_intervalbound);
					if(help_intervalbound > witness_right && witness_right_index <= intervalbounds.size()) {
						++witness_right_index;
					}
					if(help_intervalbound > witness_left + dimensional_upper_bound) {
						++witness_left_index;
					}
					LOG(INFO) << "Fall 4, addiere " << help_intervalbound;
					++j;
				}
//				LOG(INFO) << "intervalledPolynom: " << piecewisePolynoms[k-1];
//				LOG(INFO) << "tmp_intervalledPolynom: " << piecewisePolynoms[k];
				LOG(INFO) << "Witness: " << "[" << witness_left_index << ", " << witness_right_index << "]";
				LOG(INFO) << "tmp_intervalbounds: " << tmp_intervalbounds;
				LOG(INFO) << "intervalbounds: " << intervalbounds;

				std::promise<Polynom*> prom;

				piecewisePolynoms[k].push_back(tmp_intervalbounds.back(), prom.get_future() );
				jobs.emplace(k-1, witness_left_index, witness_right_index, std::move(prom));
//				DCHECK_GE(tmp_intervalledPolynom.at(tmp_intervalbounds.back())(tmp_intervalbounds.back()), 0); // Invariant: polynomial is non-negative
			}

/*#ifndef NDEBUG
			DCHECK(has_ordering(tmp_intervalbounds, std::greater<IB>())); // Invariant: the numbers of tmp_intervalbounds are strict ascendending
			for(const auto& ibound : tmp_intervalbounds) { //Invariant: the piecewise-defined polynomial is non-negative.
				DCHECK_GE(piecewisePolynoms[k].at(ibound)(ibound), 0);
			}
#endif*/
			intervalbounds.swap(tmp_intervalbounds);
			//intervalledPolynom.swap(tmp_intervalledPolynom);

			LOG(INFO) << "_old_intervals: " << intervalbounds;
		}

		constexpr size_t numthreads = 3;
		std::thread threads[numthreads];
		std::mutex queueMutex;
		auto safePop = [&] () { 
			std::lock_guard<std::mutex> guard(queueMutex);
			if(jobs.empty()) return Job(0,0,0, std::promise<Polynom*>());
			Job j(std::move(jobs.front())); 
			jobs.pop();
			return j;
		};
		auto runnable = [&] () {
			while(!jobs.empty()) { 
				Job j =std::move(safePop());
				if(j.m_witness_left_index == 0 && j.m_witness_right_index == 0 && j.m_dimension == 0) return;
				Polynom* p = sumParallelPolynomialOverWitnesses
					( dimensional_upper_bounds[j.m_dimension+1] // == dimensional_upper_bound
					  , piecewisePolynoms[j.m_dimension].bounds()
					  , SumFromZeroToUpper
					  , std::ref(piecewisePolynoms[j.m_dimension])
					  , j.m_witness_left_index
					  , j.m_witness_right_index);
				j.m_result.set_value(p);
			}};
		for(size_t i = 0; i < numthreads; ++i)
			threads[i] = std::thread(runnable);
		for(size_t i = 0; i < numthreads; ++i)
			threads[i].join();
	

		PiecedPolyAsync& intervalledPolynom = piecewisePolynoms[dimensions-1];

		LOG(INFO) << "Resulting Intervalled Polynom: " << intervalledPolynom;
#ifndef NDEBUG
		{ // checking for invariants
			//! \f$ \forall z \le \min(i_1,...,i_n) => intervalledPolynom(z) = (z+n-1 \choose z) \f$
			const unsigned int& minimalDimensionSize = *(std::min_element(dimensional_upper_bounds, dimensional_upper_bounds+dimensions));
			for(size_t z = 0; z < minimalDimensionSize; ++z) { 
				DCHECK_EQ(intervalledPolynom(z), Binomial::b(dimensions+z-1, z))  <<
					"for very small z, the upper bounds do not pose any constraint to the distribution of z." 
					"So the result is the same as for the 'bars and stars' problem.";
			}
		}
#endif
		return intervalledPolynom.extract();
	}



}//namespace

