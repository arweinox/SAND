#include "utils.hpp"

// May need to change from pass-by-reference
void RoundAngle(float& radian) {
#pragma HLS INLINE
   if (radian > PI) {
       radian -= 2 * PI;
   }
   if (radian < -PI) {
       radian += 2 * PI;
   }
}

/* Workaround for pointer reinterperation (int to float)
 * from dram bus
 *
 * dram - bus to dram from top function
 * index - location in dram for data
 * return data reinterpreted as float
 */
float ReinterpretFloat(in_32* input, int index) {
//#pragma HLS interface ap_none port=input
#pragma HLS INLINE
	din_dual reg;	// temporary register for reinterpretation
	reg.in_t = input[index];
	return reg.fp_t;
}

/*
  * Bounds the value val to the closed interval between
  * lbound and ubound
  *
  * val - value
  * lbound - lower bound
  * ubound - upper bound
  */
float BoundValue(float val, float lbound, float ubound) {
#pragma HLS INLINE region
	if (val < lbound) {
	  return lbound;
	} else if (val > ubound) {
	  return ubound;
	} else {
	  return val;
	}
}

#ifndef __SYNTHESIS__
void PrintSamples(in_32* dram) {
	for(int x = 0; x < NUM_SAMPLES; x++) {
		printf("Sample: %d\n", x);
		printf("	Object Name: %d\n", dram[RESULT_OFFSET + x*SAMPLE_SIZE + 0]);
		printf("	Density Index: %d\n", dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_INDICES(0)]);
		printf("	Pixel Index: %d\n", dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_INDICES(1)]);
		printf("	Pyramid Index: %d\n", dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_INDICES(2)]);
		printf("	Pose:  (%f, %f, %f)\n", reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJPOSE(0)]),
											reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJPOSE(1)]),
											reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJPOSE(2)]));

		printf("	Euler: (%f, %f, %f)\n", reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJEULER(0)]),
											reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJEULER(1)]),
											reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJEULER(2)]));

		printf("	Dim:   (%f, %f, %f)\n", reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJDIM(0)]),
											reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJDIM(1)]),
											reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + OBJDIM(2)]));

		printf("	Bounding box:(%d, %d), (%d, %d)\n", dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_BBOX(0)],
														dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_BBOX(1)],
														dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_BBOX(2)],
														dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_BBOX(3)]);

		printf("	Weight: %f\n\n",reinterpret_cast<float &>(dram[RESULT_OFFSET + x*SAMPLE_SIZE + SAMP_WEIGHT(3)]));
	}
}
#endif
