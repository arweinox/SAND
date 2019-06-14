#pragma once
//
// @file src/sampling/sampling.hh
// @brief base class for sampling
// @author Zhiqiang Sui
// University of Michigan, 2017
//

/* Old Headers */
//#include <renderer/scene_multi_renderer.hh>
//#include <visualization/visualization.hh>

/* New Headers */
//#include "likelihood.hpp"
//#include "inlier.hpp"
//#include "sample.hpp"
//#include "config.hpp"
//#include "LFSR.hpp"
#include "object.hpp"
//#include "utils.hpp"

namespace sampling {

/*
static double now()
{
    struct timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return tv.tv_sec + tv.tv_nsec/1e9;
}
*/

class Sampling
{
public:
    Sampling(//int num_iterations
             /*std::shared_ptr<likelihood::Likelihood> likelihood_ptr,
             std::shared_ptr<renderer::SceneMultiRenderer> ren_ptr,
             std::shared_ptr<visualization::Visualization> vis_ptr*/		 // Instead of pointers make objects in class (hopefully HLS instantiates modules)
    		);


    virtual ~Sampling() {
    	// Vivado HLS does not support dynamic calls (new, delete, malloc, alloc, and free)
    	//delete likelihood_ptr_;
    	//delete LFSR;
    }


    // void InitializePose(Sample& sample, const renderer::Heuristic& heuristic);

    // virtual void InitializePose(const renderer::ObjectName object_name, const renderer::Heuristic heuristics);

    virtual void Estimate(renderer::ObjectName obj_name) = 0;

//    virtual void get_result(Sample* top_sample) = 0;

    //virtual float get_final_weight() = 0;

    //virtual int get_final_iteration() = 0;

    void ComputeWeight();

    /* Not using this for HLS */
    //void set_likelihood_bboxes(const renderer::Heuristic heuristics[]);

protected:
    // Sample List
//  	Sample sample_list_[NUM_SAMPLES];

    // Heuristics for each sample
//    renderer::Heuristic heuristics_[NUM_SAMPLES];

  	// Weights for sample
//    Weight weights_[NUM_SAMPLES];

    // number of iterations
    //int num_iterations_;


    // Uniform RNG
//    LFSR<float> number_generator_;

    // Internal variable for weights for samples
    //Weight weight[NUM_SAMPLES];

//private:
    // Likelihood pointer
    // likelihood::Inlier likelihood_ptr_;
};


}
