/*
 * Header file contains global RAMS accessed by ILR, Inlier, and top level module.
 */
// #include "top.hpp"
#include "sandtypes.hpp"
#include "config.hpp"
#include "sample.hpp"
#include "object.hpp"
#include "objmodel_memory.hpp"
#include "obpyramid_memory.hpp"

namespace sampling {
	// Register holding current object name
	extern renderer::ObjectName curr_obj;

	// RAM holding pyramid
	extern in_32 pyramid[PYRAMID_MEM_SIZE*sizeof(in_32)];

	// RAM holding sample List for current object
  	extern Sample sample_list_[NUM_SAMPLES];

  	// RAM holding weights for each sample given a single iteration for a single object
  	extern Weight weights_[NUM_SAMPLES];
}

namespace renderer {
	// vertices of 3D model of the current object
	extern sandy::float3 object_modelv_[OBJMODEL_V_MAX];

	// faces of 3D model of the current object (holds indices of vertices in face)
	extern sandy::int3 object_modelf_[OBJMODEL_F_MAX];
}

