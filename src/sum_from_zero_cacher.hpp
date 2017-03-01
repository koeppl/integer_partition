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

namespace IntervalPartition {

    /** A simple map is used to optimize lookups of the same value. 
     * @brief Caches polynom summations \f$ \sum\limits_{k=0}^{z} p(k) \f$ for any polynom p.
     * **/
    class SumFromZeroCacher
    {
		public:
		typedef std::map<Polynom, Polynom> dict_type;

        private:
        dict_type cache; //!< caches already made queries

        public:
        /**
         * The code follows basically the proof of Lemma 4.5 with \f$ \gamma = 0 \f$.
         */
        const Polynom& operator()(const Polynom& p) {
            {
                std::map<Polynom, Polynom>::const_iterator it = cache.find(p);
                if(it != cache.end()) return it->second;
            }
            // Create a new polynomial
#ifndef NDEBUG
            auto npair = std::move(cache.emplace(p, std::move(sumFromZeroToUpper(p)) ));
            DCHECK(npair.second);
            return npair.first->second;
#else
            return cache.emplace(p, std::move(sumFromZeroToUpper(p)) ).first->second;
#endif
        }                                                                                                                                                      
    };
}//ns
#endif /* SUMFROMZEROCACHER_HPP */
