#pragma once

// #include "object.hpp"
#include "ap_int.h"


//#ifndef __SYNTHESIS__
//	typedef unsigned int in_32;
//#else
	typedef int in_32;
	#ifdef __SYNTHESIS__
		typedef ap_uint<64> in_64;
	#endif
	typedef float in_flp;	// floating point data type for input
//#endif

// Type to read DRAM bus as int or float
typedef union{
	in_32 in_t;
	float fp_t;
} din_dual;

typedef ap_uint<16> LFSR_t;

namespace sandy {

	struct int3 {
	  int a;
	  int b;
	  int c;
	};

	// GLM-like vectors
	struct vec3 {
		float v[3];
	};

	// GLM-like matrix
	struct mat4 {
		float m[4];
	};

	// CUDA-like vectors
	struct float3 {
		float x;
		float y;
		float z;
	};

	// CUDA-like vectors
	struct int4 {
		int w;
		int x;
		int y;
		int z;
	};
}
