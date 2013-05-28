/* Integer Partition
 * Computes the number of possible ordered integer partitions with upper bounds
 * Copyright (C) 2013 Roland Glück, Dominik Köppl
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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <gmpxx.h>

typedef mpz_class Z;

Z min(Z i, Z j){
    if(i<j) return i;
    return j;
}


Z height(unsigned long* bounds, int left, int right){
    Z height = 1;
    for(int i = left; i<=right; i++){
        height += bounds[i];
    }
    return height;
}

Z width(unsigned long* bounds, Z v, int left, int right){
    Z i = 0;
    Z btgwidth = 0;
    Z btgheight = height(bounds, left, right);
    if(v==0 || v==(btgheight-1)){
        return 1;
    }

    if(v<0 || v>=(btgheight)){
        return 0;
    }

    if((right-left)==0){
        return 1;
    }

    if((right-left)==1){
        Z minimum = min(bounds[left], bounds[right]);
        Z half = (bounds[left]+bounds[right])/2;
        if(v<=minimum){
            return v + 1;
        } else if(minimum < v && v <= half){
            return (min(bounds[left], bounds[right]) + 1);
        } else {
            return width(bounds, btgheight - 1 - v, left, right);
        }
    }

    if((right - left) > 1){
        int split = (right + left)/2;
        for(i=0; i<=v; i++){
            btgwidth += width(bounds, i, left, split)*width(bounds, v - i, split + 1, right);
        }
    }
    return btgwidth;
}

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
    unsigned long*const bounds = new unsigned long[bsize];
	for(size_t i = 2; i < static_cast<size_t>(argc); ++i)
		bounds[i-2] = strtoul(argv[i], NULL, 10);
	std::cout << width(bounds, z, 0, bsize-1) << std::endl;
	delete [] bounds;
	return 0;
}


