#pragma once

#include "object.hpp"
#include "weight.hpp"

// Sizes in memory (in words)
#define OBJECT_SIZE 	14
#define BBOX_SIZE 		4
#define WEIGHT_SIZE 	3
#define INDICES_SIZE 	1
#define SAMPLE_SIZE 	(OBJECT_SIZE+BBOX_SIZE+WEIGHT_SIZE+(3*INDICES_SIZE))

// Memory size for object components
#define OBJNAME_SIZE  1
#define OBJPOSE_SIZE  6
#define OBJDIM_SIZE   3
#define OBJTRANS_SIZE 4

// Macros to access Sample attributes
#define SAMP_OBJECT(X) 	X
#define SAMP_BBOX(X) 	(OBJECT_SIZE + X)	// USE XMAX,XMIN,...
#define SAMP_WEIGHT(X) 	(OBJECT_SIZE + BBOX_SIZE + X)
#define SAMP_INDICES(X) (OBJECT_SIZE + BBOX_SIZE + WEIGHT_SIZE + X)

// Access objects
#define OBJNAME 		0
#define OBJPOSE(X) 	(OBJNAME_SIZE + X)
#define OBJEULER(X)	(OBJNAME_SIZE + 3 + X)
#define OBJDIM(X)  	(OBJNAME_SIZE + OBJPOSE_SIZE + X)
#define OBJTRANS(X) (OBJNAME_SIZE + OBJPOSE_SIZE + OBJTRANS_SIZE + X)

// Access weight
#define TERMS(X) X
#define FINAL_WEIGHT(X) (X + 4)


// 24 32-bit words
struct Sample {
	// Each sample is single object
	renderer::Object object_list_;

	// Create bounding box for sample
	renderer::BoundingBox2D bbox_;

	// Weight for corresponding sample
	Weight weight_;

	//test
	float final_weight_;

	// Heatmap index (zero indexed; e.g. 0,1,2,..)
	int density_ind_;

	// Index for bounding box in heatmap indicated by density_ind_
	int pixel_ind_;

	// Flattened bouding box index, index over all heatmaps (zero indexed)
	int pyramid_heatmap_ind;
};
