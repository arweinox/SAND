//
// @file src/likelihood/balanced_outlier.cpp
// @brief Computing inliers between observation 3D points and rendered 3D points
// @author Zhiqiang Sui
// University of Michigan, 2017
//

#include "inlier.hpp"

#ifndef __SYNTHESIS__
  #include <stdio.h>
  #include <iostream>
  #include <string.h>
#endif

// computeXXX are function which use Cuda to optimize computation


using namespace likelihood;

  // Load vertices memory offsets ROM
  const int Inlier::voffset_[OBJECTS] = {
		 OBJ_V_OFFSET
  };

  // Load faces memory offsets ROM
  const int Inlier::foffset_[OBJECTS] = {
		  OBJ_F_OFFSET
  };

  // Load object vertices count ROM
  const int Inlier::vertices_num_[OBJECTS] = {
		  OBJ_V
  };

  // Load object faces count ROM
  const int Inlier::faces_num_[OBJECTS] = {
		  OBJ_F
  };


  Inlier::Inlier() : Likelihood() {}

  Inlier::~Inlier() {}

  /*
   * Calls Renderer's ResetZBuffer
   */
  void Inlier::ResetZBuffer() {
	  render_.ResetZBuffer();
  }


  /*
   * This function does the inliers between the observed scene and a single
   * hypothesized sample (object and its pose). Computation is an SSD.
   *
   * i - index of sample and weight in global RAMs
   */
   // NOTE TO SELF: Need heuristics here because of total_num_points_, otherwise
   // only need bounding_box_ from sample from computeInliersforPointCloudOnDevice
  void Inlier::ComputeInliers() //output
  {
	   for(int i = 0; i < NUM_SAMPLES; i++) {
#ifndef __SYNTHESIS__
			std::cout << "ComputeInliers.... sample " << i << " complete(Inlier)" << std::endl;
#endif


      int inliers = 0;       // accumulator, holds number of inliers
      int rendered = 0;         // count how many pixel are rendered
      int rendered_in_bbox = 0; // count how many pixel are rendered within the bounding_box

      // NOTE: Non-fixed bounds may not work during synthesis. If no, just search over
      // entire image size and only do computation if pixel coordinates fall
      // within the bounding box. In future, change so memory is size of max bounding box
      // and be smart with accessing to save memory usage
//		#ifndef __SYNTHESIS__
//			  printf("%d. Bounding box: [%d %d; %d %d]\n",index,
//					  sampling::sample_list_[index].bbox_.x_min_,
//					  sampling::sample_list_[index].bbox_.x_max_,
//					  sampling::sample_list_[index].bbox_.y_min_,
//					  sampling::sample_list_[index].bbox_.y_max_);
//		#endif

      render_.InitializeMatricies(i);
      render_.Transform(vertices_num_[sampling::curr_obj]);
      render_.Rasterize(faces_num_[sampling::curr_obj]);

		for(int x = sampling::sample_list_[i].bbox_.x_min_; x <= sampling::sample_list_[i].bbox_.x_max_; x++) {
			for(int y = sampling::sample_list_[i].bbox_.y_min_; y <= sampling::sample_list_[i].bbox_.y_max_; y++) {
			  // Grab depth of observed scene
			  float obs_point = observed_scene_pc_[WIDTH * y + x];
			  // Grab depth of pixel from renderer
			  float obj_point = render_.fetchdepth(y,x);
			  // Reset element in Z-Buffer back to FAR
			  render_.ResetElement(y,x);

			  if (obj_point > 0 && obs_point > 0) {
				  #ifdef SAME_CAMERAVIEW		// If camera view for current rendered image is same as observed scene image
					// Absolute value
					  if(obj_point -  obs_point < 0.005) { inliers++; }
				  #else
						float x_diff = 0.0f;
						float y_diff = 0.0f;
						float z_diff = 0.0f;
						// HAVE TO PUT CORECT INDICIES
						x_diff = object_pc_temp_[].x -  observed_scene_pc_[width * j + i].x;
						y_diff = object_pc_temp_[].y -  observed_scene_pc_[width * j + i].y;
						z_diff = object_pc_temp_[].z -  observed_scene_pc_[width * j + i].z;

						if (sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff) < 0.005) {    // SSD
							inliers++;
						}
				  #endif

				  rendered_in_bbox++;
				  rendered++;
				  }
			  }

//	printf("rendered: %d\n", rendered);

	// Set values for weight of sample
	#ifdef DIFFERENT_BBOXES
//		printf("Setting weights[%d];	ired:%d\n", index, rendered);
	  sampling::weights_[i].terms_[1] = rendered > 0 ? inliers / rendered : 0;		// Check if this should be rendered!!!!!!
	  sampling::weights_[i].terms_[0] = rendered_in_bbox > 0 ? inliers / rendered_in_bbox : 0;
	#else
	  // Can replace pyramid access with rendered. Does the same thing.
//	  sampling::weights_[i].terms_[0] = rendered_in_bbox > 0 ? inliers / sampling::pyramid[ATTR_CMP(FLATTEN_INDEX(sampling::sample_list_[i].density_ind_,sampling::sample_list_[index].pixel_ind_), HEURISTIC_SIZE, HEURISTICS_OFFSET, TOTAL_NUM_POINTS)] : 0;
	  sampling::weights_[i].terms_[0] = rendered_in_bbox > 0 ? inliers / rendered : 0;
	#endif
	}
		      }

}

//
//  /*
//   * Computer inlier of bounding box and observed scene (DEPRECATED)
//   */
//  void Inlier::ComputeInliersforPointCloudOnDevice()
//  {
//      // Compute the inlier for each sample
//      for(int i = 0; i < NUM_SAMPLES; i++) {
//          // Transpose object model with sample 6DOF
//          // Render Sampling for Z-Buffer
//          render_.InitializeMatricies(i);
//          render_.Transform(vertices_num_[sampling::curr_obj]);
//          render_.Rasterize(faces_num_[sampling::curr_obj]);
//
////			#ifndef __SYNTHESIS__
////				  FILE * file;
////				  file = fopen("objdepth.txt", "w");
////				  for(int x = 0; x < WIDTH; x++){
////					  for(int y =0; y < HEIGHT; y++) {
////						  fprintf(file, "%f\n", render_.fetchdepth(y,x));
////					  }
////				  }
////				 fclose(file);
////			#endif
//
//          // Compute inlier function
////          ComputeInliers(i);
//
//          // Accumulators in computeInliers do cumsum done in function below
//          //ProcessResultsForInliers(weights, num_rendered_point, num_rendered_points_in_bbox);
//      }
//  }

/*
 * Top function for computing the weights for each samples (DEPRECATED)
 *
 * depth_maps - the depth maps for each sample
 */
// void Inlier::Compute()
//  {
//      #ifdef __SYNTHESIS__
//		#ifdef VERBOSE
//        std::cout << "Computing Inliers now..." << std::endl;
//      #endif
//		#endif
//
//      #ifdef USE_NORMAL
//          #ifdef DIFFERENT_BBOXES
//              /*computeInliersforPointCloudOnDeviceWithNormal(array, observation_pc_device_, observation_normal_device_,
//                                                            horizontal_count_, vertical_count_,
//                                                            width_, height_, fx_, fy_, cx_, cy_, bounding_boxes_device_,
//                                                            weighted_arr, num_rendered_points, num_rendered_points_in_bbox);
//              */
//          #else
//              /* computeInliersforPointCloudOnDeviceWithNormal(array, observation_pc_device_, horizontal_count_, vertical_count_,
//                                                            width_, height_, fx_, fy_, cx_, cy_,
//                                                            int(bounding_box_.x_min_), int(bounding_box_.x_max_),
//                                                            int(bounding_box_.y_min_), int(bounding_box_.y_max_), weighted_arr);
//              */
//          #endif
//      #else
//          #ifdef DIFFERENT_BBOXES
////        		ComputeInliersforPointCloudOnDevice();
//          #else
//              /* computeInliersforPointCloudOnDevice(array, observation_pc_device_, horizontal_count_, vertical_count_,
//                                                  width_, height_, fx_, fy_, cx_, cy_,
//                                                  int(bounding_box_.x_min_), int(bounding_box_.x_max_),
//                                                  int(bounding_box_.y_min_), int(bounding_box_.y_max_), weighted_arr);
//              */
//          #endif
//     #endif
//  }

 /*
  * Loads BRAM with depth map of observed scene
  *
  * dram - bus to dram
  */
 void Inlier::LoadObservation(in_32 * dram) {
#pragma HLS interface m_axi port=dram bundle=master depth=99000 offset=slave

	 for(int x = 0; x < IMG_SIZE; x++) {
#pragma HLS unroll factor=2
		 observed_scene_pc_[x] = ReinterpretFloat(dram,OBSERVATION_OFFSET + x);
	 }
 }

 /*
  * Set number of faces and vertices, then load object model into global BRAM
  * of curr_obj
  *
  * dram - bus to master to load date from DDR
  */
  void Inlier::LoadObjectModel(in_32* dram) {
	  int face_offset = 0;	// Register to hold base address for face offset
	  int vert_offset = 0;  // Register to hold base address for vertices offset
 	  // Load faces
	  for(int x = 0; x < faces_num_[sampling::curr_obj]; x++){
		 face_offset = 3*x + foffset_[sampling::curr_obj];
		 renderer::object_modelf_[x].a = dram[face_offset];
		 renderer::object_modelf_[x].b = dram[face_offset+1];
		 renderer::object_modelf_[x].c = dram[face_offset + 2];
	  }
	  // Load vertices
	  for(int x = 0; x < vertices_num_[sampling::curr_obj]; x++){
 		 vert_offset = 3*x + voffset_[sampling::curr_obj];
 		 renderer::object_modelv_[x].x = ReinterpretFloat(dram, vert_offset);
 		 renderer::object_modelv_[x].y = ReinterpretFloat(dram, vert_offset + 1);
 		 renderer::object_modelv_[x].z = ReinterpretFloat(dram, vert_offset + 2);
	  }
  }
