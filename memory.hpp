#pragma once

#include "config.hpp"
#include "objmodel_memory.hpp"
#include "obpyramid_memory.hpp"

/* MEMORY OFFSETS (Data stored in this order in DDR4)***************************/
#define PYRAMID_OFFSET	   0 												// hardcoded memory address (offset from base) for observation pyramid
#define OBSERVATION_OFFSET (PYRAMID_OFFSET + PYRAMID_MEM_SIZE*NUM_OBJECTS) // hardcoded memory address (offset from base) for observations
#define OBJMODEL_OFFSET    (OBSERVATION_OFFSET + IMG_SIZE*NUM_OBJECTS) 		// hardcoded memory address (offset from base) for object point cloud
#define RESULT_OFFSET      (OBJMODEL_OFFSET + OBJMODEL_MEM_SIZE)				// hardcoded memory address (offset from base) for best sample

/* RESULTS*********************************************************************/
#define RESULTS_MEM_SIZE 7		// Memory size in words (32-bit); Pose and final_weight_
