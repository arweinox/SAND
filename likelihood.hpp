#pragma once

#include "config.hpp"
#include "sample.hpp"
#include "sandtypes.hpp"

// Problems
// 1. bounding_boxes_ and num_points_in_bounding_boxes_ are not fixed size arrays (main.json in repo sets both at 25)
// 2. xx_device_ arrays for stored on GPU

namespace likelihood {

	class Likelihood {
	public:
		Likelihood();

		virtual ~Likelihood();

		virtual void SetObservation(float* observation);

		virtual void SetObservation(sandy::float3* observation_pc);

		// Depracated
//		void SetDifferentBBoxFlag(bool different_bboxes);
//		void SetBoundingBox(const renderer::Heuristic& heuristic);
//		void SetBoundingBoxes(const renderer::Heuristic heuristics[]);
//		void SetObservationPointNumber(int point_number);
//		void SetIntrinsics(int fx, int fy, int cx, int cy);
//		void SetUseNormal(bool use_normal);
//		virtual bool Initialize() = 0;
//		virtual ~Likelihood();
//		void SetObservationNormal(sandy::float3* observation_normal);
//		void SetObservation(float* observation, int depth_channel); 	// Not used

	protected:
		// The rendered scene

/*******Deprecated************/
//		int width_;
//		int height_;
//		float fx_, fy_, cx_, cy_;
//		int vertical_count_;
//		int horizontal_count_;
//		bool different_bboxes_;
//		renderer::BoundingBox2D bounding_boxes_ [NUM_SAMPLES];

		// Might be useful
//		int depth_channel_;

		// Not used
//		int point_number_;

		// bounding box
//		renderer::Workspace bounding_box_;

		// Change to ap_uint<1> (1- true, 0 - false)
		//bool use_normal_;
		//		// Array used as actually pointer in memory for GPU (cudamalloc)
		//		sandy::float3* observation_pc_device_;
		//
		//		// Array used as actually pointer in memory for GPU (cudamalloc)
		//		sandy::float3* observation_normal_device_;

				// Array contain total_num_points_ from heuristics of each sample
				//int num_points_in_bounding_boxes_[NUM_SAMPLES];

				// observation_, observation_pc_, and observation_normal_ are set up from config file,
				// Need a pre-processing stage to set the necessary configuration parameters



		/* Pointer members */

		// Array (used for memcpy!!)
		// (size not fixed)
		float* observation_;

		// Array (used for memcpy!!)
		// (size not fixed)
		sandy::float3* observation_pc_;

		// Array (used for memcpy!!)
		// (size not fixed)
		sandy::float3* observation_normal_;

//		// Array used as actually pointer in memory for GPU (cudamalloc)
//		sandy::int4* bounding_boxes_device_;

		// Array used as actually pointer in memory for GPU (cudamalloc)
		float* observation_device_;


	};
} // namespace likelihood
