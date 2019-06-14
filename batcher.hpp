// This algorithm merges from groups of size 1 up to length of arrays
// This isn't the direct form of batcher's algorithm
// For group size p, it sorts all subgroup sizes (including p),
// starting at index j in the array, by comparing values that are k values apart

#pragma once
#include "sample.hpp"

/*
 * Swaps variable such that x < y, if y < x
 */
inline void compare(Sample &x, Sample &y) {
  if(y.weight_.final_weight_ < x.weight_.final_weight_) {
    Sample z = x;
    x = y;
    y = z;
  }
}

/*
 * Performs a indirect form of odd-even sort, which probably wont generate sorting network
 * 	- a: the array to be sorted
 * 	- r: length of the array
 */
inline void batcher(Sample a[], int n) {

	// always sorting entire array
	int l= 0;

	// length of a
	//int r = n-1;

	for (int p=1; p<n; p+=p) {    // groups of size 1,2,4,8,...
	for (int k=p; k>0; k/=2) {  // The size, of the subgroups (from size p to all halves generated from power of two); The stride, how far away from the current number to be compared
	  for (int j=k%p; j+k<n; j+=(k+k)) {  // The start position sorting of group size k
										  // Because k starts off as k=p, taking the modulus, j is equal to k (except for first time)
		for (int i=0; i<n-j-k; i++) {     // is n-j-k the size of the group (or some power of 2 division); Iterate over every element in group/subgroup
		  if ((j+i)/(p+p) == (j+i+k)/(p+p)) {  // true, if j+1 is less than a multiple of 2p and adding k does not change than answer
			compare(a[l+j+i], a[l+j+i+k]);
		  }
		}
	  }
	}
  }
}
