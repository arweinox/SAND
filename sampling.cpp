//
// @file src/sampling/sampling.cpp
// @brief base class for sampling
// @author Zhiqiang Sui
// University of Michigan, 2017
//

#include "sampling.hpp"

//#ifdef SIMULATION
// static double now() {
//      struct timespec tv;
//      clock_gettime(CLOCK_MONOTONIC, &tv);
//      return tv.tv_sec + tv.tv_nsec/1e9;
//  }
//#endif

namespace sampling
{
Sampling::Sampling(//int num_iterations
                   /*
                   std::shared_ptr<likelihood::Likelihood> likelihood_ptr,
                   std::shared_ptr<renderer::SceneMultiRenderer> ren_ptr,
                   std::shared_ptr<visualization::Visualization> vis_ptr
                   */
		){
				// Never used!
        //num_iterations_ = num_iterations;
        //iteration_ = 0;
        /*likelihood_ptr_(likelihood_ptr),
        scene_multi_renderer_ptr_(ren_ptr),
        visualization_ptr_(vis_ptr)*/
        //likelihood_ptr_;


    // Initialize likelihood
    //if (!likelihood_ptr_.Initialize())
    //{
    //	std::cerr << "Failed to initialize likelihood function" << std::endl;
    //	exit(-1);
    //}
    //std::random_device rd;
    //number_generator_;
}

/*
 * Randomly initalizes poses for NUM_SAMPLES samples
 *
 * object_name - name of current object from OBJECT objects
 * heuristics - the generic heuristics from first stage which are radnomized for each samples
 */
 // DEPRACATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*void Sampling::InitializePose(const renderer::ObjectName object_name,
							  const renderer::Heuristic heuristics)
{

	// Create all samples with random pose for given object
	for(int i = 0; i < NUM_SAMPLES; i++) {
		// Create object for sample
		renderer::Object object;

		sample_list_[i].weight_.final_weight_ = 0.0;
		sample_list_[i].bbox_ = heuristics.bbox_;
		renderer::Workspace workspace = heuristics.workspace_;

		assert(workspace.x_max_ != -std::numeric_limits<float>::infinity());
		assert(workspace.x_min_ != std::numeric_limits<float>::infinity());
		assert(workspace.y_max_ != -std::numeric_limits<float>::infinity());
		assert(workspace.y_min_ != std::numeric_limits<float>::infinity());
		assert(workspace.z_max_ != -std::numeric_limits<float>::infinity());

		// Set object name
		object.name_ = object_name;

		// Initialize pose with some random position (adhering to Gaussian distribution) (x,y,z is center of object??)
		object.pose_.pos_[0] = workspace.x_min_ + (workspace.x_max_ - workspace.x_min_) * (1.0/65535.0)*number_generator_.rand();
		object.pose_.pos_[1] = workspace.y_min_ + (workspace.y_max_ - workspace.y_min_) * (1.0/65535.0)*number_generator_.rand();
		object.pose_.pos_[2] = workspace.z_min_ + (workspace.z_max_ - workspace.z_min_) * (1.0/65535.0)*number_generator_.rand();
		object.pose_.euler_[0] = 2 * PI * (1.0/65535.0)*number_generator_.rand() - PI;
		object.pose_.euler_[1] = 2 * PI * (1.0/65535.0)*number_generator_.rand() - PI;
		object.pose_.euler_[2] = 2 * PI * (1.0/65535.0)*number_generator_.rand() - PI;

		#ifdef SIMULATION
			std::cout << object.pose_.pos_[0] << ' ' << object.pose_.pos_[1] << ' ' << object.pose_.pos_[2] << ' ' << object.pose_.euler_[0] << ' ' << object.pose_.euler_[1] << ' ' << object.pose_.euler_[2] << std::endl;
		#endif

		// I believe this is an accidennt commented out from original code
		//object.pose_ = heuristic.centroid_pose_;

		// Scaling factors
		object.dim_[0] = 1.0f;
		object.dim_[1] = 1.0f;
		object.dim_[2] = 1.0f;

		// Set randomly initialized object to sample
		sample_list_[i].object_list_ = object;
		//sample.bbox_ = heuristic.bbox_;
	}
}
*/
/*
 * Calls Inlier Compute function to compute the weights
 * for each sample.
 */
void Sampling::ComputeWeight() {
			//#ifdef SIMULATION
			//	double start = now();
			//#endif

			// Render samples on GPU and return as "array"
			// float array[NUM_SAMPLES*IMG_SIZE];
			// RenderSamples(sample_list, array);

			//#ifdef SIMULATION
			//	std::cout << "Time for rendering:                       " << (now() - start) * 1000 << " ms" << std::endl;
			//	start = now();
			//#endif

			// Call Inlier to compute inliers
			//likelihood_ptr_.Compute(sample_list_, weights_);

			//#ifdef SIMULATION
			//	std::cout << "Time for computing weight:                " << (now() - start) * 1000 << " ms" << std::endl;
			//#endif
	/*
	if (visualization_ptr_ != nullptr)
	{
			visualization_ptr_->UpdateData(array, sample_list, iteration_);
	}*/
	//scene_multi_renderer_ptr_->Unmap();
	//return weights;
}
/*
void Sampling::set_likelihood_bboxes(const renderer::Heuristic heuristics[])
{
	likelihood_ptr_.SetBoundingBoxes(heuristics);
}
*/

//void Sampling::RenderSamples(Sample sample_list[NUM_SAMPLES], float array[NUM_SAMPLES*IMG_SIZE])
//{
    /*std::vector<glm::mat4> tfs;
    renderer::TabletopScene::generateTransforms(sample_list, tfs);

    scene_multi_renderer_ptr_->UpdateTransforms(tfs);
    scene_multi_renderer_ptr_->DrawAll();
    scene_multi_renderer_ptr_->Map(&array);
		*/
//}
} // namespace sampling
