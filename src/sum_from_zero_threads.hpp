/** 
* @file sumfromzerocacher.hpp
* @brief Caches \function sumFromZeroToUpper values
* @author Dominik Köppl
* 
* @date 2015-03-06
*/

#ifndef SUMFROMZEROCACHER_HPP
#define SUMFROMZEROCACHER_HPP

#include "sum_from_zero_to_upper.hpp"
#include <shared_mutex>

namespace IntervalPartition {

    /** A simple map is used to optimize lookups of the same value. 
     * @brief Caches polynom summations \f$ \sum\limits_{k=0}^{z} p(k) \f$ for any polynom p.
     * **/
    class SumFromZeroCacherThreadSafe
    {
		public:
		typedef std::shared_timed_mutex mutex_type;
		typedef std::map<Polynom, Polynom> dict_type;

		private:
		const size_t size;
        dict_type* cache; //!< the cache caches all polynomials built with sumFromZeroToUpper
		mutex_type* mutex;

        public:
		SumFromZeroCacherThreadSafe(const size_t _size) 
			: size(_size-1)
		, cache(new dict_type[size])
		, mutex(new mutex_type[size])
		{ }
		~SumFromZeroCacherThreadSafe() {
			delete [] mutex;
			delete [] cache;
		}
        /**
         * The code follows basically the proof of Lemma 4.5 with \f$ \gamma = 0 \f$.
         */
		const Polynom& operator()(const Polynom& p) {
			DCHECK_GT(p.size(),0);
			const size_t psize = p.size()-1;
			DCHECK_LT(psize,size);
			{
				std::shared_lock<mutex_type> lock(mutex[psize]);
				dict_type::const_iterator it = cache[psize].find(p);
				if(it != cache[psize].end()) return it->second;
			}
			Polynom entry = std::move(sumFromZeroToUpper(p));
			std::lock_guard<mutex_type> lock(mutex[psize]);
			// Create a new polynom
			auto npair = std::move(cache[psize].emplace(p, std::move(entry)));
			if(npair.second) return npair.first->second; // other thread already put in some value!
			dict_type::const_iterator it = cache[psize].find(p);
			DCHECK(it != cache[psize].end());
			return it->second;
		}                                                                                                                                                      
    };
}//ns
#endif /* SUMFROMZEROCACHER_HPP */
