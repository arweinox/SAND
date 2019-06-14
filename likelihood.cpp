//
// @file src/likelihood/likelihood.cpp
// @brief Base class for computing likelihood
// @author Zhiqiang Sui
// University of Michigan, 2016
//

#include"likelihood.hpp"

//#include <cstddef>
//#include <cstring>

namespace likelihood
{
Likelihood::Likelihood() {}

Likelihood::~Likelihood() {}

void Likelihood::SetObservation(float *observation)
{
    //std::cout << "Set observation for image" << std::endl;

    observation_ = new float[IMG_SIZE];
    //std::memcpy(observation_, observation, IMG_SIZE * sizeof(float));

    for(int x = 0; x < IMG_SIZE; x++) {
		observation_[x] = observation[x];
	}
}


//void Likelihood::SetObservation(float *observation, int depth_channel)
//{
//    //std::cout << "Set observation for image" << std::endl;
//
//    /*if (observation_ != NULL)
//    {
//        delete[] observation_;;
//    }*/
//
//    /*if (observation_device_ != NULL)
//    {
//        cudaFree(observation_device_);
//    }*/
//
//    depth_channel_ = depth_channel;
//
//    observation_ = new float[IMG_SIZE * depth_channel_];
//    //std::memcpy(observation_, observation, IMG_SIZE * depth_channel_ * sizeof(float));
//
//    for(int x = 0; x < IMG_SIZE * depth_channel_; x++) {
//    	observation_[x] = observation[x];
//    }
//
//    // Copy to GPU memory
//    /*
//    cudaMalloc(&observation_device_, IMG_SIZE * depth_channel_ * sizeof(float));
//    cudaMemcpy(observation_device_, observation, IMG_SIZE * depth_channel_ * sizeof(float), cudaMemcpyHostToDevice);
//    */
//}


void Likelihood::SetObservation(sandy::float3* observation_pc)
{
    //std::cout << "Set observation for point cloud" << std::endl;

    /*if ( observation_pc_ != NULL)
    {
        delete[] observation_pc_;
    }*/

    /*if (observation_pc_device_ != NULL)
    {
        cudaFree(observation_pc_device_);
    }*/

    observation_pc_  = new sandy::float3[IMG_SIZE];
    //std::memcpy(observation_pc_, observation_pc, IMG_SIZE * sizeof(sandy::float3));

    for(int x = 0; x < IMG_SIZE; x++) {
    	observation_pc_[x] = observation_pc[x];

    }

    // Copy to GPU memory
    /*
    cudaMalloc(&observation_pc_device_, IMG_SIZE * sizeof(float3));
    cudaMemcpy(observation_pc_device_, observation_pc, IMG_SIZE * sizeof(float3), cudaMemcpyHostToDevice);
    */
}


//void Likelihood::SetObservationNormal(sandy::float3* observation_normal)
//{
//    //std::cout << "set observation for point cloud normal" << std::endl;
//
//    observation_normal_  = new sandy::float3[IMG_SIZE];
//    //std::memcpy(observation_normal_, observation_normal, IMG_SIZE * sizeof(sandy::float3));
//    for(int x = 0; x < IMG_SIZE; x++) {
//    	observation_normal_[x] = observation_normal[x];
//    }
//}

//
//void Likelihood::SetObservationPointNumber(int point_number)
//{
//    point_number_ = point_number;
//}


/*void Likelihood::SetBoundingBox(const renderer::Heuristic& heuristic)
{
    bounding_box_ = heuristic.bbox_;
}
*/

//void Likelihood::SetBoundingBoxes(const renderer::Heuristic heuristics[])
//{
//    /*if (bounding_boxes_ == nullptr)
//    {
//        bounding_boxes_ = new sandy::int4[HORIZONTAL_COUNT * VERTICAL_COUNT];
//    }
//
//    if (num_points_in_bounding_boxes_ == nullptr)
//    {
//        num_points_in_bounding_boxes_ = new int[HORIZONTAL_COUNT * VERTICAL_COUNT];
//    }
//
//    for(int i = 0; i < HORIZONTAL_COUNT * VERTICAL_COUNT; i++)
//    {
//        bounding_boxes_[i].x = int(heuristics[i].bbox_.x_min_);
//        bounding_boxes_[i].y = int(heuristics[i].bbox_.y_min_);
//        bounding_boxes_[i].z = int(heuristics[i].bbox_.x_max_);
//        bounding_boxes_[i].w = int(heuristics[i].bbox_.y_max_);
//        num_points_in_bounding_boxes_[i] = int(heuristics[i].total_num_points_);
//    }
//*/
//    // Create space and copy to GPU memory
//   /* if (bounding_boxes_device_ == nullptr)
//    {
//        cudaMalloc(&bounding_boxes_device_, HORIZONTAL_COUNT * VERTICAL_COUNT * sizeof(int4));
//    }
//
//    cudaMemcpy(bounding_boxes_device_, bounding_boxes_,
//               HORIZONTAL_COUNT * VERTICAL_COUNT * sizeof(int4), cudaMemcpyHostToDevice);*/
//}

/* Configuration file contains information */
/*
void Likelihood::SetIntrinsics(int fx, int fy, int cx, int cy)
{
    FX = fx;
    FY = fy;
    CX = cx;
    CY = cy;
}
*/

/* Configuration file contains information */
/*
void Likelihood::SetUseNormal(bool use_normal)
{
    use_normal_ = use_normal;
}
*/
} // namespace deep_scene_estimatino
