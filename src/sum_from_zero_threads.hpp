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
        private:
        std::map<Polynom, Polynom> cache; //!< caches already made queries
		std::shared_timed_mutex mutex;

        public:
        /**
         * The code follows basically the proof of Lemma 4.5 with \f$ \gamma = 0 \f$.
         */
		const Polynom& operator()(const Polynom& p) {
		{
			std::shared_lock<decltype(mutex)> lock(mutex);
			std::map<Polynom, Polynom>::const_iterator it = cache.find(p);
			if(it != cache.end()) return it->second;
		}
		std::lock_guard<decltype(mutex)> lock(mutex);
			// Create a new polynom
            auto npair = std::move(cache.emplace(p, std::move(sumFromZeroToUpper(p)) ));
			if(npair.second) return npair.first->second; // other thread already put in some value!
			std::map<Polynom, Polynom>::const_iterator it = cache.find(p);
			DCHECK(it != cache.end());
			return it->second;
        }                                                                                                                                                      
    };
}//ns
#endif /* SUMFROMZEROCACHER_HPP */
