#pragma once

#ifdef __SYNTHESIS__
	// #include "ap_int.h"
#endif
#include "global.hpp"
#include "config.hpp"
#include "memory.hpp"
#include "object.hpp"
#include "sandtypes.hpp"
#include "utils.hpp"
#include "objmodel_memory.hpp"
#include "obpyramid_memory.hpp"
#include "inlier.hpp"
#include "iterative_likelihood_reweighting.hpp"
#include <cstring>	// included for memcpy

// one-hot indicator for progress in algorithm
#define STATUS_STATES 32
#define RENDER_INIT_STATE 0x001
#define PYRAMID_STATE 0x002
#define OBJMODEL_STATE 0x004
#define INITPOSE_STATE 0x008
#define RASTERIZE_STATE 0x010
#define INLIERFUNC_STATE 0x020
#define RESAMPLE_STATE 0x040
#define DIFFUSION_STATE 0x080
#define RESET_STATE 0x100
extern ap_int<STATUS_STATES> state;


void sandtop(in_32* dram,
		ap_uint<1>* done,
		ap_uint<STATUS_STATES>* status);
