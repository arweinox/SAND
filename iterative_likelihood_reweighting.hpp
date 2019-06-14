
//
// @file src/sampling/iterative_likelihood_reweighting.hh
// @brief iterative likelihood reweighting for optimization
// @author Zhiqiang Sui
// University of Michigan, 2017
//

#pragma once

#include "sandtop.hpp"
#include "sampling.hpp"
#include "config.hpp"
#include "LFSR.hpp"
#include "utils.hpp"
#include "object.hpp"
#include "weight.hpp"
#include "batcher.hpp"
#include "sample.hpp"
#include "sandtypes.hpp"
#include "weight_helper.hpp"
#include "obpyramid_memory.hpp"

#ifndef __SYNTHESIS__
	#include <cstring>
#endif


// Marcos to load ROM
#define COEFF INLIERS_IN_BBOX, INLIERS_IN_RENDER
#define ACC_DIMS HEATMAP_1, \
			HEATMAP_1+HEATMAP_2,\
			HEATMAP_1+HEATMAP_2+HEATMAP_3,\
			HEATMAP_1+HEATMAP_2+HEATMAP_3+HEATMAP_4,\
			HEATMAP_1+HEATMAP_2+HEATMAP_3+HEATMAP_4+HEATMAP_5
#define MAP_DIM HEATMAP_H1, HEATMAP_W1,\
				HEATMAP_H2, HEATMAP_W2,\
				HEATMAP_H3, HEATMAP_W3,\
				HEATMAP_H4, HEATMAP_W4,\
				HEATMAP_H5, HEATMAP_W5


namespace sampling
{

class IterativeLikelihoodReweighting : public Sampling
{
public:
	// Instead of pointers make objects in class (hopefully HLS instantiates modules)
    IterativeLikelihoodReweighting();

    virtual void ComputeWeight();

    virtual bool CheckConverge();

    virtual bool CheckRestart();

    virtual void Diffusion();

    void DiffuseTranslation(int index, int coord, float trans_std);

    void TransferSamples(Sample samples1[NUM_SAMPLES], Sample samples2[NUM_SAMPLES], int index1, int index2);


    void CalculateAccWeightsFromDensityPyramid();

    void Estimate(renderer::ObjectName obj_name);

    // Changed from vector of heuristics to arrays
    virtual void InitializePose(int count, renderer::ObjectName obj_name);

    virtual void Resample();

    virtual void RestartSample();

    void RestartSampleFromPreviousRun();

    void IncrementIteration();

    float get_final_weight();

    // int get_final_iteration();

	// iteration counter
    int iteration_;
    // restart counter
    int restart_times_;

    // Depracated (HLS wont synthesis. Internal variables set in Estimate.
//    void InitEstimate(renderer::ObjectName obj_name);

    // Returns samples with highest weight for object
//    void get_result(Sample* top_sample);

    /********************************** For testing only */
    void get_acc_weight(float * weight_bram);
    void set_weights(float weight[NUM_SAMPLES]);
    void get_weights(float * weights);
    void set_samples(Sample * samples[NUM_SAMPLES]);
    void get_sample_list(Sample * samples_out);

    /*********************************************************/

    // Should add function to preform cum sum, used in Resample and InitializePose

protected:

    Sample most_likely_sample_for_each_restart_[NUM_SAMPLES];

    // inliers_in_box is the first index, inlier_in_render is the second
    static const float coefficients_[2];

    // Only need a single normal distribution since standard deviation changes every loop
    static const float normal_dist[1000];

    // Random number generator for InitializePose
    LFSR<float> number_generator_;

    // Accumulator for cumsum on the number of pixels on each level of pyramid
    static const int acc_dims_[MAX_PYRAMID_HEATMAPS];

    // Accumulator for cumsum of pixel densities(weights)
    float acc_weights_[PYRAMID_PIXELS];

    // Object name of the current object
    renderer::ObjectName obj_name_;

	// Map dimensions (32 bit width follows height)
	static const int map_dim_[10];
};
}
