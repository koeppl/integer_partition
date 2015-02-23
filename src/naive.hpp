/** 
 * @file naive.hpp
 * @brief Naive computation algorithm for ordered integer partitions
 * @author Dominik Köppl
 * 
 * @date 2015-02-23
 */
#ifndef NAIVE_HPP
#define NAIVE_HPP

#include <algorithm>
#include <cassert>
#include <glog/logging.h>

/** 
 * Calculates the integer partition of z with upper bounds[left..right]
 * 
 * @param bounds The maximum contraints i_1,...,i_n
 * @param v the integer that shall be partitioned
 * @param left left interval index of bounds
 * @param right right interval index of bounds. Note that right is *inclusive*.
 * @pre \f$ 0 \le left \le right \f$ must hold
 * @pre \f$ right < bounds.length \f$ must hold
 * 
 * @return the number of integer partitions
 */
template<class t_Z>
t_Z naive_bounds(const unsigned int*const bounds, t_Z v, size_t left, size_t right){
    t_Z i = 0;
    t_Z summedBound = std::accumulate(bounds+left, bounds+right+1, 0);

    if(v==0 || v==summedBound) return 1;
    if(v<0 || v>summedBound) return 0;
    if(right == left) return 1;

    if(right == 1+left) {
        const t_Z minimum = std::min(bounds[left], bounds[right]);
        const t_Z arith = (bounds[left]+bounds[right])/2;
		if(v > arith) v = summedBound-v;
		return std::min(v,minimum)+1;
    }

	DCHECK_GT(right, 1+left);
    t_Z ret = 0;
	const size_t split = (right + left)/2;
	for(i=0; i<=v; i++){
		ret += naive_bounds<t_Z>(bounds, i, left, split)*naive_bounds<t_Z>(bounds, v - i, split + 1, right);
	}
	DCHECK_GT(ret, 0);
    return ret;
}

#endif /* NAIVE_HPP */
