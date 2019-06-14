#include "sandtop.hpp"

// Define global RAMs
namespace sampling {
	// Register holding current object name
	renderer::ObjectName curr_obj;

	// RAM holding pyramid
	in_32 pyramid[PYRAMID_MEM_SIZE*sizeof(in_32)];

	// RAM holding sample List for current object
	Sample sample_list_[NUM_SAMPLES];

	// RAM holding weights for each sample given a single iteration for a single object
	Weight weights_[NUM_SAMPLES];
}
namespace renderer {
	// vertices of 3D model of the current object
	sandy::float3 object_modelv_[OBJMODEL_V_MAX];

	// faces of 3D model of the current object (holds indices of vertices in face)
	sandy::int3 object_modelf_[OBJMODEL_F_MAX];
}


/*
 * Performs incomplete iteration and pulls the sample with the greatest weight
 *
 * top_sample - sample with highest weight
 * dram - bus to external RAM
 */
void sandtop(in_32* dram, ap_uint<1>* done, ap_uint<STATUS_STATES>* status) {
#pragma HLS interface m_axi port=dram bundle=master depth=99000 offset=slave

#pragma HLS interface s_axilite port=done  bundle=slave offset=0x010
#pragma HLS interface s_axilite port=status bundle=slave offset=0x020
#pragma HLS interface s_axilite port=return bundle=slave

	// ROM holding the names for each particle
	const renderer::ObjectName object_list[OBJECTS+1] = {
			renderer::blue_cup,
			renderer::clorox,
			renderer::coke,
			renderer::detergent,
			renderer::downy,
			renderer::ranch,
			renderer::red_bowl,
			renderer::salt,
			renderer::scotch_brite,
			renderer::spray_bottle,
			renderer::sugar,
			renderer::sunscreen,
			renderer::tide,
			renderer::toy,
			renderer::waterpot,
			renderer::none
	};

	// Instantiate likelihood and sampling
	likelihood::Inlier likelihood;
	sampling::IterativeLikelihoodReweighting sampler;


	// Register to hold memory address for results of current object
	int results_offset = 0;
	// Register to hold memory address for observation pyramid of current object
	int pyramid_offset = PYRAMID_OFFSET;

	// Load observation scene into BRAM
	likelihood.LoadObservation(dram);
	// Reset entire Z-Buffer
	likelihood.ResetZBuffer();



	for(int obj = 0; obj < NUM_OBJECTS; obj++) {

		// Load observation pyramid into global BRAM
		memcpy(sampling::pyramid, dram, PYRAMID_MEM_SIZE*sizeof(in_32));

		// Set current object
		sampling::curr_obj = object_list[obj];

		// Load object model
		likelihood.LoadObjectModel(dram);

		// Initialize bounding box and workspace of samples
		sampler.InitializePose(NUM_SAMPLES, object_list[obj]);
//	sampler.InitializePose(NUM_SAMPLES, object_list[0]);


		#ifndef __SYNTHESIS__
			printf("Initial Samples>>>>>>>>>>>>\n");
			memcpy( dram + RESULT_OFFSET, (in_32 *)sampling::sample_list_, NUM_SAMPLES*SAMPLE_SIZE*sizeof(float));
			PrintSamples(dram);
			std::cout << "Entering loop: " << std::endl;
		#endif

//		while (!sampler.CheckConverge() && sampler.iteration_ < NUM_ITERATIONS) // Iteration conditional my be redundant
//		{
			// Inlier rasterize object model
			likelihood.ComputeInliers();

			// Compute Likelihood function to assess accuracy of guess
			sampler.ComputeWeight();

			// Resample to get more samples with better accuracy
			sampler.Resample();
#ifndef __SYNTHESIS__

			printf("After Resample\n\n");
			memcpy( dram + RESULT_OFFSET, (in_32 *)sampling::sample_list_, NUM_SAMPLES*SAMPLE_SIZE*sizeof(float));
			PrintSamples(dram);
#endif


			// Pertub samples 6DOF to get variety amongst samples
			sampler.Diffusion();


//			sampler.iteration_ += 1;
//
//			#ifndef __SYNTHESIS__
//				printf("Iteration %d finished\n\n", sampler.iteration_);
//				std::cout << "Iteration: " << sampler.iteration_ + 1 << std::endl;
//
//			memcpy( dram + RESULT_OFFSET, (in_32 *)sampling::sample_list_, NUM_SAMPLES*SAMPLE_SIZE*sizeof(float));
//			PrintSamples(dram);
//			#endif
//		}
//
//		// Inlier Compute Weight function call
//		likelihood.ComputeInliers();
//		// Compute weight for the last round
//		sampler.ComputeWeight();
//
//		#ifndef __SYNTHESIS__
//			std::cout << "Number of iterations: " << sampler.iteration_ << std::endl;
//		#endif
//
//		// Calculate address to write results
//		results_offset = RESULT_OFFSET + RESULTS_MEM_SIZE * object_list[obj];

//	    for(int x=0; x < NUM_SAMPLES; x++) {
//	    	dram[results_offset] = sampling::sample_list_[x].weight_.final_weight_;
//	    	dram[results_offset+1] = sampling::sample_list_[x].object_list_.pose_.pos_[0];
//	    	dram[results_offset+2] = sampling::sample_list_[x].object_list_.pose_.pos_[1];
//	    	dram[results_offset+3] = sampling::sample_list_[x].object_list_.pose_.pos_[2];
//	    	dram[results_offset+4] = sampling::sample_list_[x].object_list_.pose_.euler_[0];
//			dram[results_offset+5] = sampling::sample_list_[x].object_list_.pose_.euler_[1];
//			dram[results_offset+6] = sampling::sample_list_[x].object_list_.pose_.euler_[2];
//	    }

//	    // Increment observation pyramid offset
//	    pyramid_offset += PYRAMID_MEM_SIZE;
//
//	    in_32* results = dram + RESULT_OFFSET;
//
//		//For now send entire sample list (only need 6DOF and final_weight_)
//		memcpy(results, (in_32 *)sampling::sample_list_, NUM_SAMPLES*SAMPLE_SIZE*sizeof(float));
//
	}
//
//	// Indicate module is done
//	done[0] = (ap_uint<1>)1;
//	#ifndef __SYNTHESIS__
//		std::cout << "Done!" << std::endl;
//	#endif
}
