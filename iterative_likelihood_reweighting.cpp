//
// @file src/sampling/iterative_likelihood_reweighting.cpp
// @brief iterative likelihood reweighting for optimization
// @author Zhiqiang Sui
// University of Michigan, 2017
//

//#define CONFIGURU_IMPLEMENTATION 1

#include "iterative_likelihood_reweighting.hpp"

// All values for pass to function in software are fixed class member variables
// Normal distributions should be calculate before running, stored in array and randomly accessed
// Uniform distributed value are calculated from LFSR and scaled appropriately for the specific case
// Parameters for Likelihood pointer are hardcoded!!!!!!!!!

// Current problem: reassigning values in the member arrays (sample_list and heuristics_)

using namespace sampling;

// Load cumsum on heatmap dimensions ROM
const int IterativeLikelihoodReweighting::acc_dims_[MAX_PYRAMID_HEATMAPS] = {
		ACC_DIMS
};
// Load coefficients for final weight equation ROM
const float IterativeLikelihoodReweighting::coefficients_[2] = {
		COEFF
};
// Load map dimensions ROM
const int IterativeLikelihoodReweighting::map_dim_[10] = {
		MAP_DIM
};
// Load normal distribution ROM
const float IterativeLikelihoodReweighting::normal_dist[1000] = {
	#include "normal_dist.dat"
};

/*
 * Constructor sets the default object to none, restart counter to 0, and initializes LFSR
 */
IterativeLikelihoodReweighting::IterativeLikelihoodReweighting() : obj_name_(renderer::coke), restart_times_(0) // Need to fix this
{
	number_generator_.set_seed(6);		// Seed LFSR
}

/*
 * Driving function for particle filter and main function for ILR. Performs
 * weight computation, resampling, diffusion to find
 */
void IterativeLikelihoodReweighting::Estimate(renderer::ObjectName obj_name)
{
    // Set iteration_ to 0
    iteration_ = 0;
    restart_times_ = 0;
    obj_name_ = obj_name;

    // Initialize bounding box and workspace of samples
    InitializePose(NUM_SAMPLES, obj_name);

    while (CheckConverge() == false || iteration_ < NUM_ITERATIONS)
    {
		#ifndef __SYNTHESIS__
        	std::cout << "Iteration: " << iteration_ + 1 << std::endl;
		#endif

        // Compute Liklihood function to assess accuracy of guess
        ComputeWeight();
        // Resample to get more samples with better accuracy
        Resample();
        // Pertub samples 6DOF to get variety amongst samples
        Diffusion();

        iteration_ += 1;
    }

    // Compute weight for the last round
    ComputeWeight();

	#ifndef __SYNTHESIS__
    	std::cout << "Number of iterations: " << iteration_ << std::endl;
  //  std::cin.get();
	#endif
}


/*
 * Forms normalized cumulative sum as an array of densities from the first stage
 */
void IterativeLikelihoodReweighting::CalculateAccWeightsFromDensityPyramid()
{
    // flattened observation density stack
    float sum_density;

    sum_density = pyramid[ATTR(0, DENSITY_SIZE, DENSITIES_OFFSET)];
    acc_weights_[0] = pyramid[ATTR(0, DENSITY_SIZE, DENSITIES_OFFSET)];
    for(int c = 1; c < PYRAMID_PIXELS; c++) {
        sum_density += pyramid[ATTR(c, DENSITY_SIZE, DENSITIES_OFFSET)];
        acc_weights_[c] = pyramid[ATTR(c, DENSITY_SIZE, DENSITIES_OFFSET)] + acc_weights_[c-1];
    }

    // Normalize densities
    for(int norm_ind = 0; norm_ind < PYRAMID_PIXELS; norm_ind++) {
        acc_weights_[norm_ind] /= sum_density;
    }

	#ifndef __SYNTHESIS__
	  std::cout << "	weight cum sum done" << std::endl;
	#endif
}

/*
 * Randomly initialize NUM_SAMPLES samples across entire observation denisty pyramid
 *
 * count - number of samples to initialize
 * obj_name - name of object
 */
void IterativeLikelihoodReweighting::InitializePose(int count, renderer::ObjectName obj_name)
{
	#ifndef __SYNTHESIS__
		std::cout << ("=========================Initialize pose...") << std::endl;
	#endif
	// Set iteration_ to 0
	iteration_ = 0;
	restart_times_ = 0;
	obj_name_ = obj_name;

    CalculateAccWeightsFromDensityPyramid();

    for(int sample_ind = 0; sample_ind < count; sample_ind++)
    {
        int density_ind = -1;   // heat map index
        int pixel_ind = -1;     // pixel index
        int sampled_flattened_index = 0;
        float rand_number = number_generator_.uniform_dist();		// Get random number from uniform distribution

        // Randomly select pixel from pyramid
        for(int i = 0; i < PYRAMID_PIXELS; i++) {
            if (rand_number < acc_weights_[i])
            {
                sampled_flattened_index = i;
                i =  PYRAMID_PIXELS;   // This stops the loop; This may be a problem
            }
        }

        // Calculate density index and the corresponding pixel index
        for(int i = 0; i < PYRAMID_HEATMAPS; i++)
        {
            if (sampled_flattened_index < acc_dims_[i])  // May be problem because you can't unroll
            {
                density_ind = i;
                if ( i > 0 )
                {
                    pixel_ind = sampled_flattened_index - acc_dims_[i-1];
                }
                else
                {
                    pixel_ind = sampled_flattened_index;
                }

                break;
            }
        }


        //assert(density_ind < PYRAMID_HEATMAPS);


        //if (density_pyramid_[density_ind].heuristics_[pixel_ind].total_num_points_ < 100)
        if (pyramid[ATTR_CMP(sampled_flattened_index, HEURISTIC_SIZE, HEURISTICS_OFFSET, TOTAL_NUM_POINTS)] < 100)
        {
#ifndef __SYNTHESIS__
			std::cout << "Initialize Pose number points miss: " << std::endl;
#endif
            sample_ind--;
            //continue;
        } else {

        	//printf("Initialize sample\n");
          //assert(pixel_ind >= 0);
          //assert(pixel_ind < acc_dims_[density_ind]-(density_ind==0 ? 0 : acc_dims_[density_ind-1]));

			#ifndef __SYNTHESIS__
//				printf("sample index: %d\n", sampled_flattened_index);
			#endif

			// Create object for sample
			renderer::Object object;
			//renderer::Workspace workspace = pyramid_.heuristics[sampled_flattened_index].workspace_;
			int workspace_offset = ATTR_CMP(sampled_flattened_index, HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(XMIN));
			int bbox_offset = ATTR_CMP(sampled_flattened_index, HEURISTIC_SIZE, HEURISTICS_OFFSET, BBOX(XMIN));

//#ifndef __SYNTHESIS__
//	        std::cout << "density_ind: " << density_ind <<
//	        		"	pixel ind: " << pixel_ind <<
//					"	bbox_offset: " << pyramid[bbox_offset] << std::endl;
//#endif


			object.name_ = obj_name_;
			sample_list_[sample_ind].final_weight_ = -99.0;
			sample_list_[sample_ind].density_ind_ = density_ind;
			sample_list_[sample_ind].pixel_ind_ = pixel_ind;
			sample_list_[sample_ind].pyramid_heatmap_ind = sampled_flattened_index;
			sample_list_[sample_ind].bbox_.x_min_ = pyramid[bbox_offset];
			sample_list_[sample_ind].bbox_.x_max_ = pyramid[++bbox_offset];
			sample_list_[sample_ind].bbox_.y_min_ = pyramid[++bbox_offset];
			sample_list_[sample_ind].bbox_.y_max_ = pyramid[++bbox_offset];

			// Replace call with contents of function
			//Sampling::InitializePose(sample_ind, heuristics[sampled_flattened_index]);	// Behavior unsure, since sample_list_ and heuristics_ are the member variables of base class Sampling

			//assert(workspace.x_max_ != -std::numeric_limits<float>::infinity());
			//assert(workspace.x_min_ != std::numeric_limits<float>::infinity());
			//assert(workspace.y_max_ != -std::numeric_limits<float>::infinity());
			//assert(workspace.y_min_ != std::numeric_limits<float>::infinity());
			//assert(workspace.z_max_ != -std::numeric_limits<float>::infinity());

			// Extract floats as workspace
			float xmin_temp = ReinterpretFloat(pyramid, workspace_offset + XMIN);
			float xmax_temp = ReinterpretFloat(pyramid, workspace_offset + XMAX);
			float ymin_temp = ReinterpretFloat(pyramid, workspace_offset + YMIN);
			float ymax_temp = ReinterpretFloat(pyramid, workspace_offset + YMAX);
			float zmin_temp = ReinterpretFloat(pyramid, workspace_offset + ZMIN);
			float zmax_temp = ReinterpretFloat(pyramid, workspace_offset + ZMAX);

			// Initialize pose with some random position (adhering to Gaussian distribution) (x,y,z is center of object??)
			object.pose_.pos_[0] = xmin_temp + (xmax_temp - xmin_temp) * number_generator_.uniform_dist();
			object.pose_.pos_[1] = ymin_temp + (ymax_temp - ymin_temp) * number_generator_.uniform_dist();
			object.pose_.pos_[2] = zmin_temp + (zmax_temp - zmin_temp) * number_generator_.uniform_dist();
			object.pose_.euler_[0] = 2 * PI * number_generator_.uniform_dist() - PI;
			object.pose_.euler_[1] = 2 * PI * number_generator_.uniform_dist() - PI;
			object.pose_.euler_[2] = 2 * PI * number_generator_.uniform_dist() - PI;

//			#ifndef __SYNTHESIS__
//			int two = 1;
//			int three = 3;
//				float rnum = number_generator_.uniform_dist();
//				printf("random num: %f\n", rnum);
//				printf("xmin: %x\n", pyramid[workspace_offset + XMIN]);
//				printf("xmax: %x\n", pyramid[workspace_offset + XMAX]);
//				printf("difference: %d\n", pyramid[workspace_offset + XMAX] - pyramid[workspace_offset + XMIN]);
//				printf("workspace xmin: %f\n",
//						(int)(pyramid[workspace_offset + XMAX] - pyramid[workspace_offset + XMIN])*rnum);
//			#endif

			// I believe this is an accident commented out from original code
			//object.pose_ = heuristic.centroid_pose_;

			// Scaling factors
			object.dim_[0] = 1.0f;
			object.dim_[1] = 1.0f;
			object.dim_[2] = 1.0f;

			// Set randomly initialized object to sample
			sample_list_[sample_ind].object_list_ = object;
		}
    }
}

/*
 * Computes the weight for all samples through Inlier,
 * calculates the final weight, then sort the samples
 * by final weight.
 */
void IterativeLikelihoodReweighting::ComputeWeight()
{
//#pragma HLS inline region
//#pragma HLS interface s_axilite port=weights register register_mode=forward

#ifndef __SYNTHESIS__
	std::cout << "=========================Compute Weight..." << std::endl;
#endif


		// Call Inlier to compute inliers
//		Inlier::Compute(sample_list_, weights_);

		// Compute final weight for each sample; weight from all objects of sample contribute to final weight
		for(size_t i = 0; i < NUM_SAMPLES; i++) {
			float total_weight = 0.0;			// temporary variable
			Weight weight = weights_[i];
			for (int j = 0; j < 2; j++) {		// UNROLL
				//#pragma HLS UNROLL factor=2
				#ifndef __SYNTHESIS__
					assert(weight.terms_[j] >= 0);
					assert(weight.terms_[j] <= 1);
				#endif

				total_weight += coefficients_[j] * weight.terms_[j];

				#ifdef ENABLE_SIGMOID
					//assert(Sigmoid(-NORM_RANGE_FOR_SIGMOID + weight.terms_[j] * 2 * NORM_RANGE_FOR_SIGMOID) >= 0);
					//assert(Sigmoid(-NORM_RANGE_FOR_SIGMOID + weight.terms_[j] * 2 * NORM_RANGE_FOR_SIGMOID) <= 1);
					total_weight += coefficients_[j] * Sigmoid(-NORM_RANGE_FOR_SIGMOID + weight.terms_[j] * 2 * NORM_RANGE_FOR_SIGMOID);
				#endif
			}

			#ifdef ENABLE_SIGMOID
				total_weight = Sigmoid(-NORM_RANGE_FOR_SIGMOID + total_weight * 2 * NORM_RANGE_FOR_SIGMOID);
			#endif

			// Set weight parameters for sample
			weight.final_weight_ = total_weight;
			sample_list_[i].weight_ = weight;
		}

		// If Batcher Sorting Network does not work out, store a NUM_SAMPLES array
		// where the index corresponds to the sample and contains a struct of the pose_estimation
		// and the weight (sort on weights)
    // Sort
    batcher(sample_list_, NUM_SAMPLES);

    //std::cout << "Maximum weight: " << sample_list_.back().weight_.final_weight_ << std::endl;
    //std::cout << "Maximum weight (inliner in render): " << sample_list_.back().weight_.terms_["inliers_in_render"] << std::endl;
}

///
// Convergence is limited by (1) restart times and (2) iterations
bool IterativeLikelihoodReweighting::CheckConverge()
{
//#pragma HLS inline region
    if (CONVERGENCE_THRESHOLD > 0.0)
    {
        if (restart_times_ > MAX_RESTART_TIMES)
        {
            return true;
        }

        return sample_list_[NUM_SAMPLES - 1].weight_.final_weight_ >= CONVERGENCE_THRESHOLD;
    }

    return iteration_ >= NUM_ITERATIONS;
}

bool IterativeLikelihoodReweighting::CheckRestart()
{
#pragma HLS INLINE

    // Restart orientation of the sample each restart_per_iteration iteration
    if ((iteration_ + 1) % RESTART_PER_ITERATION == 0)
    {
        most_likely_sample_for_each_restart_[restart_times_] = sample_list_[NUM_SAMPLES - 1];
        restart_times_ += 1;

//        if (restart_times_ > MAX_RESTART_TIMES)
//        {
//            return false;
//        }


        return !(restart_times_ > MAX_RESTART_TIMES);
    }
    return false;
}


// Randomaly shifting, over Gaussian distribution. the particles dedicated to a given bounding box slightly to give unique possibilities for resampling
void IterativeLikelihoodReweighting::Diffusion()
{
//#pragma HLS INLINE region

#ifndef __SYNTHESIS__
	std::cout << "=========================Diffusion..." << std::endl;
#endif

    int iteration_remainder = iteration_ % RESTART_PER_ITERATION;

    // should combine into one if-statement
   	//if (!(CONVERGENCE_THRESHOLD > 0 && !CheckRestart())) // CheckRestart returns false is restart_times_ > MAX_RESTART_TIMES == true; second statement is redundant
	//{

	// Normal distributions are now computed in pre-stage and stored
	float decay = 1 + iteration_remainder * iteration_remainder * DECAY;

	// Normal distribution generated every call
	float trans_std =  TRANSITION_VARIANCE / decay;
	float rot_std =  ROTATION_VARIANCE / decay;

	//int sample_ind = 0;
	for(int i = 0; i < NUM_SAMPLES - 1; i++)
	{
//			//Sample sample = sample_list_[i];
//			//renderer::Workspace workspace = heuristics_[i].workspace_;
//			int workspace_offset = ATTR_CMP(FLATTEN_INDEX(sample_list_[i].density_ind_,sample_list_[i].pixel_ind_), HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(XMIN));
//			// Sample contains single object, no need to for-each loop
//			renderer::Object object = sample_list_[i].object_list_;

		// Generate random number from uniform distribution [0, 12]
		float rand_op = number_generator_.rand() * 13/65535;

#ifndef __SYNTHESIS__
			printf("random_operation: %f\n", rand_op);
#endif

// Switch conditional to int 0 - 12
		if ( rand_op < 1 )
		{
			DiffuseTranslation(i, POSE_X, trans_std);
		}
		else if ( rand_op < 2 )
		{
			DiffuseTranslation(i, POSE_Y, trans_std);
		}
		else if ( rand_op < 3 )
		{
			sample_list_[i].object_list_.pose_.euler_[2] += normal_dist[number_generator_.index_normal()]/rot_std;

			RoundAngle(sample_list_[i].object_list_.pose_.euler_[2]);
		}
		else if ( rand_op < 4 )
		{
			DiffuseTranslation(i, POSE_Z, trans_std);
		}
		else if ( rand_op < 5 )
		{
			sample_list_[i].object_list_.pose_.euler_[0] += normal_dist[number_generator_.index_normal()]/rot_std;

		}
		else if ( rand_op < 6 )
		{
			sample_list_[i].object_list_.pose_.euler_[1] += normal_dist[number_generator_.index_normal()]/rot_std;

			RoundAngle(sample_list_[i].object_list_.pose_.euler_[1]);
		}
		else if ( rand_op < 7 )
		{
			DiffuseTranslation(i, POSE_X, trans_std);
			DiffuseTranslation(i, POSE_Y, trans_std);
			DiffuseTranslation(i, POSE_Z, trans_std);
		}
		else if ( rand_op < 8 )
		{
			sample_list_[i].object_list_.pose_.euler_[0] += normal_dist[number_generator_.index_normal()]/rot_std;
			sample_list_[i].object_list_.pose_.euler_[1] += normal_dist[number_generator_.index_normal()]/rot_std;
			sample_list_[i].object_list_.pose_.euler_[2] += normal_dist[number_generator_.index_normal()]/rot_std;

			RoundAngle(sample_list_[i].object_list_.pose_.euler_[0]);
			RoundAngle(sample_list_[i].object_list_.pose_.euler_[1]);
			RoundAngle(sample_list_[i].object_list_.pose_.euler_[2]);
		}
		else if ( rand_op < 9 )
		{
			sample_list_[i].object_list_.pose_.euler_[0] += PI;
			RoundAngle(sample_list_[i].object_list_.pose_.euler_[0]);
		}
		else if ( rand_op < 10 )
		{
			sample_list_[i].object_list_.pose_.euler_[1] += PI;
			RoundAngle(sample_list_[i].object_list_.pose_.euler_[1]);
		}
		else if ( rand_op < 11 )
		{
			sample_list_[i].object_list_.pose_.euler_[2] += PI;
			RoundAngle(sample_list_[i].object_list_.pose_.euler_[2]);
		}
		else if (rand_op < 12)
		{
			#ifdef SIMULATION
				std::cout << "Diffusing among current density map" << std::endl;
			#endif

//					int map_width = pyramid_.map_dim[sample.density_ind_].width;
//					int map_height = pyramid_.map_dim[sample.density_ind_].height;
//					int map_width = pyramid[ATTR_CMP(sample_list_[i].density_ind_, MAP_DIM_SIZE, MAP_DIM_OFFSET, MAP_WIDTH)];
//					int map_height = pyramid[ATTR_CMP(sample_list_[i].density_ind_, MAP_DIM_SIZE, MAP_DIM_OFFSET, MAP_HEIGHT)];

			// Get map dimensions
			int map_index = 2*sample_list_[i].density_ind_;
			int map_width = map_dim_[map_index];
			int map_height = map_dim_[map_index+1];

			int x_ind = sample_list_[i].pixel_ind_ % map_width;
			int y_ind = sample_list_[i].pixel_ind_ / map_width;  // is this right?

			// Randomly decide to change either the x or y coordinate +/- 1
			float direction = number_generator_.uniform_dist();
			if (direction <= 0.25)
			{
				y_ind = y_ind - 1 >= 0 ? y_ind - 1 : 0;
			}
			else if (direction <= 0.5)
			{
				y_ind = y_ind + 1 < map_height ? y_ind + 1 : map_height - 1;
			}
			else if (direction <= 0.75)
			{
				x_ind = x_ind - 1 >= 0 ? x_ind - 1 : 0;
			}
			else if (direction <= 1 )
			{
				x_ind = x_ind + 1 < map_width ? x_ind + 1 : map_width - 1;
			}

			#ifndef __SYNTHESIS__
//				assert(y_ind >= 0);
//				assert(y_ind < map_height);
//				assert(x_ind >= 0);
//				assert(x_ind < map_width);
			#endif

			int pixel_ind_temp =  y_ind * map_width + x_ind;

			// Compute offsets for access pyramid
			int bbox_offset = ATTR_CMP(sample_list_[i].pyramid_heatmap_ind, HEURISTIC_SIZE, HEURISTICS_OFFSET, BBOX(XMIN));
			int workspace_offset = ATTR_CMP(sample_list_[i].pyramid_heatmap_ind, HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(XMIN));

#ifndef __SYNTHESIS__
			printf("bbox_offset: %d workspace_offset: %d", bbox_offset, workspace_offset);
#endif

			// Note: density map is 1D array, pixel_ind in the index in the 1D array for sample
			sample_list_[i].pixel_ind_ = pixel_ind_temp;
			sample_list_[i].bbox_.x_min_ = pyramid[bbox_offset];
			sample_list_[i].bbox_.x_max_ = pyramid[++bbox_offset];
			sample_list_[i].bbox_.y_min_ = pyramid[++bbox_offset];
			sample_list_[i].bbox_.y_max_ = pyramid[++bbox_offset];

			// If any coordinate exceeds dimension than set to max
			if (sample_list_[i].object_list_.pose_.pos_[0] < ReinterpretFloat(pyramid, workspace_offset))
			{
					sample_list_[i].object_list_.pose_.pos_[0] = ReinterpretFloat(pyramid, workspace_offset);
			}
			else if (sample_list_[i].object_list_.pose_.pos_[0] > ReinterpretFloat(pyramid, workspace_offset + XMAX))
			{
					sample_list_[i].object_list_.pose_.pos_[0] = ReinterpretFloat(pyramid, workspace_offset + XMAX);
			}

			if (sample_list_[i].object_list_.pose_.pos_[1] < ReinterpretFloat(pyramid, workspace_offset + YMIN))
			{
					sample_list_[i].object_list_.pose_.pos_[1] = ReinterpretFloat(pyramid, workspace_offset + YMIN);
			}
			else if (sample_list_[i].object_list_.pose_.pos_[1] > ReinterpretFloat(pyramid, workspace_offset + YMAX))
			{
					sample_list_[i].object_list_.pose_.pos_[1] = ReinterpretFloat(pyramid, workspace_offset + YMAX);
			}

			if (sample_list_[i].object_list_.pose_.pos_[2] < ReinterpretFloat(pyramid, workspace_offset + ZMIN))
			{
					sample_list_[i].object_list_.pose_.pos_[2] = ReinterpretFloat(pyramid, workspace_offset + ZMIN);
			}
			else if (sample_list_[i].object_list_.pose_.pos_[2] > ReinterpretFloat(pyramid, workspace_offset + ZMAX))
			{
					sample_list_[i].object_list_.pose_.pos_[2] = ReinterpretFloat(pyramid, workspace_offset + ZMAX);
			}
		}
		// Diffuse across density maps
		else if (rand_op < 13)
		{
			std::cout << "Diffusing across density maps" << std::endl;

			int new_density_ind = 0;
			int density_ind = sample_list_[i].density_ind_;	// Temporary register for quick access

			// Get map dimensions
			int map_index = 2*sample_list_[i].density_ind_;
			int map_width =  map_dim_[map_index];
			int map_height = map_dim_[map_index+1];

			// Randomly select previous or next heatmap
			float direction = number_generator_.uniform_dist();
			if (direction <= 0.5){
				new_density_ind = density_ind - 1 >= 0 ? density_ind - 1 : PYRAMID_HEATMAPS - 1;
			}
			else if (direction <= 1) {
				new_density_ind = density_ind + 1 < PYRAMID_HEATMAPS ? density_ind + 1 : 0;
			}


				#ifndef __SYNTHESIS__
//				assert(new_density_ind >= 0);
//				assert(new_density_ind < PYRAMID_HEATMAPS);
				#endif

			// Get the pixel index with aspect ratio
			int x_ind = sample_list_[i].pixel_ind_ % map_width;
			int y_ind = sample_list_[i].pixel_ind_ / map_width;
			// Get new pixel index
			int new_y_ind = y_ind * map_dim_[2*new_density_ind + MAP_HEIGHT] / map_height;
			int new_x_ind = x_ind * map_dim_[2*new_density_ind + MAP_WIDTH] / map_width;

#ifndef __SYNTHESIS__

			std::cout << "i: " << i << "  pixel_ind:_ " << sample_list_[i].pixel_ind_  << std::endl;
			std::cout << "x_ind: " << x_ind;
			std::cout << " 	y_ind: " << y_ind;
			std::cout << "  new_y_ind: " << new_y_ind;
			std::cout << " 	new_x_ind: " << new_x_ind << std::endl;
#endif

#ifndef __SYNTHESIS__
//					assert(new_y_ind >= 0);
//					assert(new_y_ind < pyramid_.map_dim[new_density_ind].height);
//					assert(new_x_ind >= 0);
//					assert(new_x_ind < pyramid_.map_dim[new_density_ind].width);
#endif


			// Register to hold pixel_ind to prevent continuous re-reads
			int pixel_temp = new_y_ind * map_dim_[2*new_density_ind + MAP_WIDTH] + new_x_ind;
			sample_list_[i].pixel_ind_ = pixel_temp;

			// Get new density index
			sample_list_[i].density_ind_ = new_density_ind;

			// Flatten bounding box index over all heatmaps
			int flattened_ind = sample_list_[i].pyramid_heatmap_ind;
#ifndef __SYNTHESIS__
std::cout << "new_density_ind: " << new_density_ind <<
		"	pixel_temp: " << pixel_temp  <<
		"		acc_dims:  " <<
		acc_dims_[new_density_ind-1] << std::endl;
#endif
			int bbox_offset = ATTR_CMP(flattened_ind, HEURISTIC_SIZE, HEURISTICS_OFFSET, BBOX(XMIN));
			int workspace_offset = ATTR_CMP(flattened_ind, HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(XMIN));

#ifndef __SYNTHESIS__
			std::cout << "flattened_ind: " <<  flattened_ind << "	bbox_offset: " << bbox_offset << "	workspace_offset: " << workspace_offset << std::endl;
#endif

			sample_list_[i].bbox_.x_min_ = pyramid[bbox_offset];
			sample_list_[i].bbox_.x_max_ = pyramid[++bbox_offset];
			sample_list_[i].bbox_.y_min_ = pyramid[++bbox_offset];
			sample_list_[i].bbox_.y_max_ = pyramid[++bbox_offset];


//					renderer::Workspace workspace =  heuristics_[acc_dims_[new_density_ind] + sample_list_[i].pixel_ind_].workspace_;

			// If any coordinate exceeds dimension than set to max
			if (sample_list_[i].object_list_.pose_.pos_[0] < ReinterpretFloat(pyramid, workspace_offset))
			{
					sample_list_[i].object_list_.pose_.pos_[0] = ReinterpretFloat(pyramid, workspace_offset);
			}
			else if (sample_list_[i].object_list_.pose_.pos_[0] > ReinterpretFloat(pyramid, workspace_offset + XMAX))
			{
					sample_list_[i].object_list_.pose_.pos_[0] = ReinterpretFloat(pyramid, workspace_offset + XMAX);
			}

			if (sample_list_[i].object_list_.pose_.pos_[1] < ReinterpretFloat(pyramid, workspace_offset + YMIN))
			{
					sample_list_[i].object_list_.pose_.pos_[1] = ReinterpretFloat(pyramid, workspace_offset + YMIN);
			}
			else if (sample_list_[i].object_list_.pose_.pos_[1] > ReinterpretFloat(pyramid, workspace_offset + YMAX))
			{
					sample_list_[i].object_list_.pose_.pos_[1] = ReinterpretFloat(pyramid, workspace_offset + YMAX);
			}

			if (sample_list_[i].object_list_.pose_.pos_[2] < ReinterpretFloat(pyramid, workspace_offset + ZMIN))
			{
					sample_list_[i].object_list_.pose_.pos_[2] = ReinterpretFloat(pyramid, workspace_offset + ZMIN);
			}
			else if (sample_list_[i].object_list_.pose_.pos_[2] > ReinterpretFloat(pyramid, workspace_offset + ZMAX))
			{
					sample_list_[i].object_list_.pose_.pos_[2] = ReinterpretFloat(pyramid, workspace_offset + ZMAX);
			}
		}
	}
}

/*
 * Helper function for Diffuse that perturbs position for translation
 *
 * index - index of sample
 * coord - indicator of which coordinate to perturb (x=0, y=1, z=2)
 * trans_std - standard deviation for random number used for perturbation
 */
void IterativeLikelihoodReweighting::DiffuseTranslation(int index, int coord, float trans_std)
{
//#pragma HLS inline region

	float disturb;
	float curr_val = sample_list_[index].object_list_.pose_.pos_[coord];		// current coordinate value
//    float disturb = normal_dist_trans[(int)(number_generator_.rand() * 1000/65535)]/trans_std;
    int offset = ATTR_CMP(sample_list_[index].pyramid_heatmap_ind, HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(2*coord));
    float range_min = ReinterpretFloat(pyramid, offset);
	float range_max = ReinterpretFloat(pyramid, offset + 1);

    int repeat_ind = 0;		// loop counter

    // Restrict inside the workspace
    while ((curr_val + disturb < range_min || curr_val + disturb > range_max) && repeat_ind < 20)
    {
        disturb = normal_dist[number_generator_.index_normal()]/trans_std;

#ifndef __SYNTHESIS__
			printf("index: %d   normal distribution: %f\n",number_generator_.index_normal(), normal_dist[number_generator_.index_normal()]);
#endif

        repeat_ind++;
    }

    if (repeat_ind == 20)
    {
        if (curr_val + disturb < range_min)
        {
        	sample_list_[index].object_list_.pose_.pos_[coord] = range_min;
        }
        if (curr_val + disturb > range_max)
        {
        	sample_list_[index].object_list_.pose_.pos_[coord] = range_max;
        }
    }
    else
    {
    	sample_list_[index].object_list_.pose_.pos_[coord] += disturb;
    }
}

// Randomly sample over existing samples over Gaussian distribution.
// Called every Estimate call
void IterativeLikelihoodReweighting::Resample()
{
//#pragma HLS inline region

#ifndef __SYNTHESIS__
	std::cout << "=========================Resample..." << std::endl;
#endif

    // Normalization
    float weight_sum = 0.0f;
    for(int x = 0; x < NUM_SAMPLES; x++) {
    	// Unroll
        weight_sum += sample_list_[x].weight_.final_weight_;
    }

    float normalized_weights[NUM_SAMPLES];
    for(int x = 0; x < NUM_SAMPLES; x++) {
    	// Unroll
        normalized_weights[x] = (sample_list_[x].weight_.final_weight_ / weight_sum);
    }

    // Cumulative sum of normalized weights vector
    float acc_weights[NUM_SAMPLES];
    acc_weights[1] = normalized_weights[0];
    for(int x = 1; x < NUM_SAMPLES; x++) {
    	acc_weights[x] = normalized_weights[x] + acc_weights[x-1];
    }

    Sample resampled_sample_list[NUM_SAMPLES];
    int num_sample = 0;
    while (num_sample < NUM_SAMPLES-1)    // Randomly sample over existing samples (but for the samples that have a lower probability?? You're selected from sample_list_ for the first few samples which correspond to sample with lower weigths)
    {						// I think it's replace the first n number of samples with lowest probability


        // Uniform randomly generate number [0, 1]
        float rand_number = number_generator_.uniform_dist();

        int sample_index = 0;
        for(size_t sample_ind = 0; sample_ind < NUM_SAMPLES; sample_ind++)
        {
            if (rand_number < acc_weights[sample_ind])
            {
                sample_index = sample_ind;
                break;
            }
        }
        TransferSamples(sample_list_,resampled_sample_list, sample_index, num_sample);
        //resampled_sample_list[num_sample] = sample_list_[sample_index];
        num_sample++;
    }

    //Keep the most likely particle which is the last one in the list
//    resampled_sample_list[NUM_SAMPLES - 1] = sample_list_[NUM_SAMPLES - 1];
    TransferSamples(sample_list_,resampled_sample_list, NUM_SAMPLES - 1, NUM_SAMPLES - 1);

    // Store
    for(int x = 0; x < NUM_SAMPLES; x++) {
         //sample_list_[x] = resampled_sample_list[x];
         TransferSamples(resampled_sample_list, sample_list_, x, x);
    }
}

/*
 * Transfer samples attributes from samples1 to sample2
 *
 * samples1 - array of Samples
 * samples2 - array of Samples
 * index1 - index to access samples1
 * index2 - index to access samples2
 */
void IterativeLikelihoodReweighting::TransferSamples(Sample* samples1, Sample* samples2, int index1, int index2){
//#pragma HLs interface ap_bus port=samples1
//#pragma HLs interface ap_bus port=samples2
#pragma HLS INLINE
	samples2[index2].object_list_.name_ = samples1[index1].object_list_.name_;

	samples2[index2].object_list_.dim_[0] = samples1[index1].object_list_.dim_[0];
	samples2[index2].object_list_.dim_[1] = samples1[index1].object_list_.dim_[1];
	samples2[index2].object_list_.dim_[2] = samples1[index1].object_list_.dim_[2];

	samples2[index2].object_list_.pose_.pos_[0] = samples1[index1].object_list_.pose_.pos_[0];
	samples2[index2].object_list_.pose_.pos_[1] = samples1[index1].object_list_.pose_.pos_[1];
	samples2[index2].object_list_.pose_.pos_[2] = samples1[index1].object_list_.pose_.pos_[2];

	samples2[index2].object_list_.pose_.euler_[0] = samples1[index1].object_list_.pose_.euler_[0];
	samples2[index2].object_list_.pose_.euler_[1] = samples1[index1].object_list_.pose_.euler_[1];
	samples2[index2].object_list_.pose_.euler_[2] = samples1[index1].object_list_.pose_.euler_[2];


	samples2[index2].bbox_.x_max_ = samples1[index1].bbox_.x_max_;
	samples2[index2].bbox_.x_min_ = samples1[index1].bbox_.x_min_;
	samples2[index2].bbox_.y_max_ = samples1[index1].bbox_.y_max_;
	samples2[index2].bbox_.y_min_ = samples1[index1].bbox_.y_min_;

	samples2[index2].weight_ = samples1[index1].weight_;
	samples2[index2].density_ind_ = samples1[index1].density_ind_;
	samples2[index2].pixel_ind_ = samples1[index1].pixel_ind_;
	samples2[index2].pyramid_heatmap_ind = samples1[index1].pyramid_heatmap_ind;
}


/*
 * Randomly reinitialize pose and scaling of the first NUM_RESTART_SAMPLE samples.
 *  (At restart either re-initialize each sample or update sample list
 *  by changing the first num_restart_sample (the first num_restart_sample unlikely [reminder: sample_list is sorted]))
 *  samples with a perturbed copy of the most likely sample over all restarts
 *  Called from 'CheckRestart' which is called in 'Diffusion' call in 'Estimate'
 */
void IterativeLikelihoodReweighting::RestartSample()
{
//#pragma HLS inline region
    // If restart times reach the max_restart_times, re-initialize the samples from the
    // most likely sample from the previous run
//    if (restart_times_ == MAX_RESTART_TIMES)
//    {
//        RestartSampleFromPreviousRun();
//
//    } else {
		// Reinitialize NUM_RESTART_SAMPLE samples 6DOF pose and scaling
		for(int i = 0; i < NUM_RESTART_SAMPLE; i++) {

				// Replication of code in InitializePose (DO SOMETHING ABOUT THIS)
				renderer::Object object;
	//			renderer::Workspace workspace = pyramid_.heuristics[i].workspace_;

				int workspace_offset = ATTR_CMP(i, HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(XMIN));
	//			int bbox_offset = ATTR_CMP(i, HEURISTIC_SIZE, HEURISTICS_OFFSET, BBOX(XMIN));

	//			assert(workspace.x_max_ != -std::numeric_limits<float>::infinity());
	//			assert(workspace.x_min_ != std::numeric_limits<float>::infinity());
	//			assert(workspace.y_max_ != -std::numeric_limits<float>::infinity());
	//			assert(workspace.y_min_ != std::numeric_limits<float>::infinity());
	//			assert(workspace.z_max_ != -std::numeric_limits<float>::infinity());

				// Set object name
				object.name_ = obj_name_;

				// Extract floats as workspace
				float xmin_temp = ReinterpretFloat(pyramid, workspace_offset + XMIN);
				float xmax_temp = ReinterpretFloat(pyramid, workspace_offset + XMAX);
				float ymin_temp = ReinterpretFloat(pyramid, workspace_offset + YMIN);
				float ymax_temp = ReinterpretFloat(pyramid, workspace_offset + YMAX);
				float zmin_temp = ReinterpretFloat(pyramid, workspace_offset + ZMIN);
				float zmax_temp = ReinterpretFloat(pyramid, workspace_offset + ZMAX);

				// Initialize pose with some random position (adhering to Gaussian distribution) (x,y,z is center of object??)
				object.pose_.pos_[0] = xmin_temp + (xmax_temp - xmin_temp) * number_generator_.uniform_dist();
				object.pose_.pos_[1] = ymin_temp + (ymax_temp - ymin_temp) * number_generator_.uniform_dist();
				object.pose_.pos_[2] = zmin_temp + (zmax_temp - zmin_temp) * number_generator_.uniform_dist();
				object.pose_.euler_[0] = 2 * PI * number_generator_.uniform_dist() - PI;
				object.pose_.euler_[1] = 2 * PI * number_generator_.uniform_dist() - PI;
				object.pose_.euler_[2] = 2 * PI * number_generator_.uniform_dist() - PI;

	//			#ifndef __SYNTHESIS__
	//				std::cout << object.pose_  .pos_.x << ' ' << object.pose_.pos_.y << ' ' << object.pose_.pos_.z << ' ' << object.pose_.euler_.x << ' ' << object.pose_.euler_.y << ' ' << object.pose_.euler_.z << std::endl;
	//			#endif


				// I believe this is an accident commented out from original code
				//object.pose_ = heuristic.centroid_pose_;

				// Scaling factors
				object.dim_[0] = 1.0f;
				object.dim_[1] = 1.0f;
				object.dim_[2] = 1.0f;

				// Set randomly initialized object to sample
				sample_list_[i].object_list_ = object;
//		}
    }
}

/*
 *  Uses current sample_list for next restart. Updates list by replacing the
 *  first NUM_RESTART_SAMPLE (the first NUM_RESTART_SAMPLE unlikely [
 *  reminder: sample_list is sorted]) samples with a perturbed copy of the
 *  most likely sample over all restarts
 */
void IterativeLikelihoodReweighting::RestartSampleFromPreviousRun()
{

	#ifndef __SYNTHESIS__
		std::cout << "restart from the most likely sample from the previous run" << std::endl;
	#endif

    //Sort
    batcher(most_likely_sample_for_each_restart_, NUM_SAMPLES);

    for(int i = 0; i < NUM_RESTART_SAMPLE; i++)
    {
        // Set for num_restart_samples to the mostly_likely samples from all restarts
    	sample_list_[i] = most_likely_sample_for_each_restart_[NUM_SAMPLES - 1];

		sample_list_[i].object_list_.pose_.pos_[0] += number_generator_.uniform_dist();
		sample_list_[i].object_list_.pose_.pos_[1] += number_generator_.uniform_dist();
		sample_list_[i].object_list_.pose_.pos_[2] += number_generator_.uniform_dist();

		sample_list_[i].object_list_.pose_.euler_[0] += number_generator_.uniform_dist();
		sample_list_[i].object_list_.pose_.euler_[1] += number_generator_.uniform_dist();
		sample_list_[i].object_list_.pose_.euler_[2] += number_generator_.uniform_dist();
    }
}

/*
 * Return weight of most likely sample
 */
float IterativeLikelihoodReweighting::get_final_weight()
{
//#pragma HLS inline region
    return sample_list_[NUM_SAMPLES - 1].weight_.final_weight_;
}

/* TESTING******************************************************************/
///*
// *  For testing Resample() and RestartSample()
// */
//void IterativeLikelihoodReweighting::set_weights(float weight[NUM_SAMPLES]) {
//#pragma HLS interface ap_bus port=weight depth=4938
//	for(int x = 0; x < NUM_SAMPLES; x++){
//		sample_list_[x].weight_.final_weight_  = weight[x];
//	}
//}

///*
// *  For testing Resample() and RestartSample()
// */
//void IterativeLikelihoodReweighting::set_samples(Sample * samples[NUM_SAMPLES]) {
//	for(int x = 0; x < NUM_SAMPLES; x++){
//		sample_list_[x] = samples[x];
//	}
//}

///*
// *  For testing Resample() and RestartSample()
// */
//void IterativeLikelihoodReweighting::get_acc_weight(float * weight_bram) {
//	#pragma HLS interface ap_bus port=weight_bram depth=4938
//	for(int x = 0; x < PYRAMID_PIXELS; x++) {
//		weight_bram[x] = acc_weights_[x];
//	}
//}
//
//void IterativeLikelihoodReweighting::get_sample_list(Sample * samples_out) {
//#pragma HLS interface ap_bus port=samples_out depth=25000
//	for(int x = 0; x < NUM_SAMPLES; x++) {
//		*(samples_out + x) = sample_list_[x];
//	}
//}
//
//void IterativeLikelihoodReweighting::get_weights(float * weights){
//#pragma HLS interface ap_bus port=weights depth=625
//	for(int x = 0; x < NUM_SAMPLES; x++) {
//		*(weights + x) = sample_list_[x].weight_.final_weight_;
//	}
//}

/***********************DEPRACATED OR USELESS CODE **************************************/
/****************************************************************************************/

// Increments the iteration count
/*void IterativeLikelihoodReweighting::IncrementIteration() {
	iteration_++;
}*/



/*int IterativeLikelihoodReweighting::get_final_iteration()
{
//#pragma HLS inline region
    return iteration_;
} */

/*
 * Initialize sampler internal variables with the object name and the
 * observation density pyramid along with tracking variables (HLS is
 * not going to synthesize this)
 *
 * obj_name - name of the current object
 */
//void IterativeLikelihoodReweighting::InitEstimate(renderer::ObjectName obj_name) {
//	iteration_ = 0;
//	restart_times_ = 0;
//	obj_name_ = obj_name;
//}
//// Returns result of rewerghting max of either possible samples from current iteration
//// or for all time
//void IterativeLikelihoodReweighting::get_result(Sample* top_sample)
//{
//#pragma HLS interface ap_bus port=top_sample
////#pragma HLS inline region
//    if (CONVERGENCE_THRESHOLD > 0)
//    {
//
//        /*Sample max_sample_iter =  std::max_element(most_likely_sample_for_each_restart_.begin(),
//                                most_likely_sample_for_each_restart_.end(),
//                                [] (const Sample& a, const Sample& b)  {
//                                    return a.weight_.final_weight_ < b.weight_.final_weight_;
//                                });*/
//    	// Find sample with max final weight
//    	Sample max_sample_iter = most_likely_sample_for_each_restart_[0];
//        for (int x = 0; x < NUM_SAMPLES; x++) {
//        	if(most_likely_sample_for_each_restart_[x].weight_.final_weight_ > max_sample_iter.weight_.final_weight_) {
//        		max_sample_iter = most_likely_sample_for_each_restart_[x];
//        	}
//        }
//
//        /*for(const auto& sample : most_likely_sample_for_each_restart_)
//        {
//            std::cout << sample.weight_.final_weight_ << std::endl;
//        }*/
//
//        if (max_sample_iter.weight_.final_weight_ > sample_list_[NUM_SAMPLES - 1].weight_.final_weight_)
//        {
//        	*top_sample = max_sample_iter;
//        }
//    }
//
//    top_sample = (sample_list_ + NUM_SAMPLES - 1);
//}
/*
 * Initalizes the pose for NUM_SAMPLES samples of object. All samples initialized
 * with same heuristic before InitEstimate() call.
 *
 * object_list - a single object belonging to a sample
 * heuristics - heuristics given as input for object type from first stage
 */
/*
void IterativeLikelihoodReweighting::InitializePose(const renderer::ObjectName object_name,
                                                    const renderer::Heuristic heuristics)
{
    //sample_list_.clear();

    // Initialize particles
    //for(int i = 0; i < NUM_SAMPLES; i++)
    //{
        //Sample sample;
        //sample.weight_.final_weight_ = 0.0;
        //sample.object_list_  = object_list;
        //sample.bbox_ = heuristics.bbox_;

        Sampling::InitializePose(object_name, heuristics);

        // Check that the position exist in the bounds of image/scene/workspace
        assert(sample.object_list_.pose_.pos_[0] >= heuristics.workspace_.x_min_);
        assert(sample.object_list_.pose_.pos_[0] <= heuristics.workspace_.x_max_);

        assert(sample.object_list_.pose_.pos_[1] >= heuristics.workspace_.y_min_);
        assert(sample.object_list_.pose_.pos_[1] <= heuristics.workspace_.y_max_);

        assert(sample.object_list_.pose_.pos_[2] >= heuristics.workspace_.z_min_);
        assert(sample.object_list_.pose_.pos_[2] <= heuristics.workspace_.z_max_);

        //sample_list_[i] = sample;
        //heuristics_[i] = heuristics[i];
  // }

    //std::memcpy(heuristics_, heuristics, sizeof(renderer::Heuristic));
}
*/
