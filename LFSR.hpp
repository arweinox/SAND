/* Simple 16-bit LFSR based off implementation on Wikipedia */
#pragma once

#include "ap_int.h"
#include "sandtypes.hpp"

template<class T>
class LFSR {
public:

	LFSR() : seed_(0) {}

	LFSR(LFSR_t seed) : seed_(seed) {}

	// Generate a random 16-bit number
	T rand() {
//	 #pragma HLS interface ap_bus port=num
		static LFSR_t state = this->seed_;
		//uint32_t start_state = 0xACE1u;  /* Any nonzero start state will work. */

		LFSR_t  lfsr = state;
		LFSR_t  bit;                    /* Must be 16bit to allow bit<<15 later in the code */
		LFSR_t  period = 0;
		LFSR_t  nbits = 16;					// LFSR length (data types only support 16 bits)

		for(int x = 0; x < 16; x++){
			/* taps: 16 14 13 11; feedback polynomial: x^16 + x^14 + x^13 + x^11 + 1 */
			bit  = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 5) ) & 1;
			lfsr =  (lfsr >> 1) | (bit << 15);
		}

		// Set state at current output for next call
		state = lfsr;
		return (T)lfsr;
	}

	// Generate random number uniformly distributed 0.0 to 1.0
	T uniform_dist() {
		return (1.0/65535.0)*rand();
	}

	// Generate random number uniformly distributed min to max
	T uniform_dist(float min, float max) {
		return (max-min)*rand()/65535.0 + min;
	}

	// Index normal distribution
	int index_normal() {
		return (1000/65535)*rand();
	}

	// Set the seed value
	void set_seed(LFSR_t seed) {
		this->seed_ = seed;
	}

private:
		LFSR_t seed_;

	// Believe was for testing, used from site which created the code
	/*float uniform_dist() {
		float min = 1.0;
		float max = 4294967295.0;

		return (static_cast<float>(rand_()) - min)/(max - min);
	}*/
};
