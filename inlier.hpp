#pragma once
//
// @file src/likelihood/inlier.hh
// @brief Computing inliers between observation 3D points and rendered 3D points
// @author Zhiqiang Sui
// University of Michigan, 2017
//

#include "object.hpp"
#include "global.hpp"
#include "memory.hpp"
#include "config.hpp"
#include "renderer.hpp"
#include "likelihood.hpp"
#include "obpyramid_memory.hpp"
#include "objmodel_memory.hpp"


#define OBJ_V  OBJ_1V, OBJ_2V, OBJ_3V, OBJ_4V, OBJ_5V, \
			   OBJ_6V, OBJ_7V, OBJ_8V, OBJ_9V, OBJ_10V, \
			   OBJ_11V, OBJ_12V, OBJ_13V, OBJ_14V, OBJ_15V

#define OBJ_F  OBJ_1F, OBJ_2F, OBJ_3F, OBJ_4F, OBJ_5F, \
			   OBJ_6F, OBJ_7F, OBJ_8F, OBJ_9F, OBJ_10F, \
			   OBJ_11F, OBJ_12F, OBJ_13F, OBJ_14F, OBJ_15F

#define OBJ_V_OFFSET OBJ_1V_OFFSET, OBJ_2V_OFFSET, OBJ_3V_OFFSET, \
						OBJ_4V_OFFSET, OBJ_5V_OFFSET, OBJ_6V_OFFSET, OBJ_7V_OFFSET, \
						OBJ_8V_OFFSET, OBJ_9V_OFFSET, OBJ_10V_OFFSET, OBJ_11V_OFFSET, \
						OBJ_12V_OFFSET, OBJ_13V_OFFSET, OBJ_14V_OFFSET, OBJ_15V_OFFSET

#define OBJ_F_OFFSET OBJ_1F_OFFSET, OBJ_2F_OFFSET, OBJ_3F_OFFSET, OBJ_4F_OFFSET, OBJ_5F_OFFSET, OBJ_6F_OFFSET, OBJ_7F_OFFSET, OBJ_8F_OFFSET, OBJ_9F_OFFSET, OBJ_10F_OFFSET, OBJ_11F_OFFSET, OBJ_12F_OFFSET, OBJ_13F_OFFSET, OBJ_14F_OFFSET, OBJ_15F_OFFSET

namespace likelihood {

class Inlier : public Likelihood
{
public:

		Inlier();

		~Inlier();

		void ResetZBuffer();

//		// Macros if-statement. Only calls ComputeInliersforPointCloudOnDevice
//		virtual void Compute();

		// Compute inlier for a single sample
		virtual void ComputeInliers();

		// Iterate through samples and renders the object with 6DOF of sample, calls ComputesInliers
//		virtual void ComputeInliersforPointCloudOnDevice();

		// Load the depth map of the observed scene
		virtual void LoadObservation(in_32 * dram);

		// Load object model faces and vertices into global BRAM
		virtual void LoadObjectModel(in_32* dram);

private:
	// Depth map of the observed scene
	float observed_scene_pc_[IMG_SIZE];

	// Z-buffer stored in render and inlier accesses it
//	// Z-buffer for the current sample
//	float object_z_buffer_[MAX_OBJECT_PC];

	// Renderer for transformation and z-buffer generation
	renderer::Renderer<float> render_;

	// Number of faces for current object
	static const int faces_num_[OBJECTS];

	// Number of vertices for current object
	static const int vertices_num_[OBJECTS];

	// Memory offset for each object's vertices
	static const int voffset_[OBJECTS];

	// Memory offset for each object's face
	static const int foffset_[OBJECTS];
};
} // namespace likelihood
