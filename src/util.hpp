#ifndef UTIL_HPP
#define UTIL_HPP

/** 
 * @file util.hpp
 * @brief Utility Functions
 * @author Dominik Köppl
 * 
 * @date 2015-03-03
 */


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

}//ns

#include <ostream>
#include <vector>
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	os << "(";
	for(size_t i = 0; i < v.size(); ++i) {
		os << v[i] << " ";
	}
	os << ")";
	return os;
}

#include <functional>
#include <map>
/* From http://slackito.com/2011/03/17/automatic-memoization-in-cplusplus0x/ */

/*
template <typename ReturnType, typename... Args>
std::function<ReturnType (Args...)> memoize(std::function<ReturnType (Args...)> func)
{
    return ([=](Args... args) mutable {
    		static std::map<std::tuple<Args...>, ReturnType> cache;
            std::tuple<Args...> t(args...);
            if (cache.find(t) == cache.end())                
                cache[t] = func(args...);
            return cache[t];
    });
}
*/





template <typename R, typename... Args>
std::function<const R& (Args...)> memoize(R (*fn)(Args...)) {
	std::map<std::tuple<Args...>, R> cache;
	return [fn,cache](Args... args) mutable -> const R& {

		auto argt = std::make_tuple(args...);
		{
			auto memoized = cache.find(argt);
			if(memoized != cache.end()) return memoized->second;
		}
		auto npair = std::move(cache.emplace(argt, std::move(fn(args...))));
		DCHECK(npair.second);
		return npair.first->second; 
	};
}

#endif /* UTIL_HPP */
