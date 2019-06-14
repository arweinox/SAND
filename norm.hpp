#include "LFSR.hpp"


/*
 * Computes a normal distribution using CLT by adding 12 LFSR
 * generated random sequences together, then selecting a random
 * value from the collectoin with the normal distribution using
 * a LFSR generated value.
 */

uint16_t norm(int mu, int std) {
    LFSR<float> urng;      // Uniform LFSR random number generator
    int len = 100;  // len of the random numbers for each sequence
    unsigned normal[len] = {0}; // Accumulator for summation of random number sequences

    // Random seeds
    uint16_t seeds[] = {0x002, 0xACE1, 0x0B0D,
                        0x234, 0x0042, 0x7AB1,
                        0x9A0, 0x0892, 0x049E,
                        0xF2A, 0x0693, 0xC892};

    for(int i = 0; i < 12; i++){
      urng.set_seed(seeds[i]);
      for(int x = 0; x <len; x++) {
        normal[x] = normal[x] + urng.rand();
      }
    }

    // Calculate current mu (since values are large mean > len)
    int mean = 0;
    for(int x = 0; x < len; x++) {
      mean = mean + normal[x];
    }
    mean = mean/len;

    int std_est = 0;
    // Calculate standard error
    for(int x = 0; x < len; x++) {
      std_est = std_est + pow(normal[x] - mean, 2);
    }
    std_est = std_est / (len-1);

    // Change distribution for input parameters
    for(int x = 0; x < len; x++) {
      normal[x] = (normal[x] - mean + mu)*std/std_est;
    }

    // NEED TO CHANGE THIS SO THAT IT STORES DISTRIBUTION IN BRAM

    // return random value from array with normal distribution
    //urng.set_seed(0xbeef);
    return normal[urng.rand()];
}
