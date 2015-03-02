
//! to pass complex elements with commas as parameters to a macro, add wrap the parameter with this macro
#define MACRO_ESCAPE(...) __VA_ARGS__ 

#define TO_STR 	(x)  #x
#define TO_XSTR (x)  TO_STR(x)

/** Concatenates by preprocessor directive **/
#define CONCATENATE_DETAIL(x, y) x##y
/** Concatenates by preprocessor directive **/
#define CONCATENATE(x, y) CONCATENATE_DETAIL(x, y)
/** Creates an unique identifier for a given string **/
#define MAKE_UNIQUE(x) CONCATENATE(x, __COUNTER__)

/** Returns the bounds of any container **/
#define BOUNDS(x) x.begin(), x.end()

/** Internal use **/
#define forI_(counterVar, counterS, boundS, func) \
	{ \
		for(size_t counterS = 0; counterS < boundS; ++counterS) { \
			const size_t& counterVar = counterS; \
			{ \
				func \
			} \
		} \
	}

/** 
 * Generates a \code for(size_t counterVar = 0; counterVar < bound; ++counterVar) { func } \endcode loop
 * 
 * @param counterVar the counter-variable used in func
 * @param bound the upper bound of counterVar
 * @param func the codeblock in the for-loop
 */
#define forI(counterVar, bound, func) forI_(counterVar, MAKE_UNIQUE(jVar), bound, func)

/** Internal use **/
#define forIC_(counterVar, boundVar, bound, func) \
	{ \
		const size_t boundVar = bound; \
		forI_(counterVar, MAKE_UNIQUE(j), boundVar, func) \
	}

/** 
 * Generates a \code const size_t size = bound; for(size_t counterVar = 0; counterVar < size; ++counterVar) { func } \endcode loop
 * 
 * @param counterVar the counter-variable used in func
 * @param bound the upper bound of counterVar, stored in a variable.
 * @param func the codeblock in the for-loop
 */
#define forIC(counterVar, bound, func)  forIC_(counterVar, MAKE_UNIQUE(bVar), bound, func)

