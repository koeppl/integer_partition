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


/*
#include <cassert>
#include <vector>
#include <map>
#include <iostream>
#include <cstdint>
#include <stdint.h>
#include <cmath>
#include <algorithm>

*/
#include "interval_partition.hpp"


// output handling
/*
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
	template<class U, class V>
	std::ostream& operator<<(std::ostream& os, const std::pair<U, V>& pair)
	{
		os << "|";
		os << pair.first << " -> " << pair.second;
		os << "|";
		return os;
	}

*/





int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << argv[0] << " - calculate the " << std::endl;
		std::cout << "Usage: " << argv[0] << " n i_0 [i_1 [i_2 [...]]]" << std::endl;
		return 1;
	}
	const size_t bsize = argc-2;
	const unsigned long z = strtoul(argv[1], NULL, 10);
    unsigned int*const bounds = new unsigned int[bsize];
	for(size_t i = 2; i < static_cast<size_t>(argc); ++i)
		bounds[i-2] = strtoul(argv[i], NULL, 10);

	IntervalPartition::IntervalledPolynom intervalledPolynom = IntervalPartition::generateIntervalPartition(bounds, bsize);
	std::cout << intervalledPolynom(z) << std::endl;
	delete [] bounds;
	return 0;
}

//	long sum = 0;
//	for(size_t i = 0; i < bsize;++i) sum += bounds[i];
//	std::cout << intervalledPolynom(sum/2) << std::endl;



