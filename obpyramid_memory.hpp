#pragma once
/*
 * This header files contains macros for accessing
 * and manipulating observation pyramid stored
 * within the DDR4.
 */

//=======================================/
//
//	Macros for addressing pyramid
//
//  Ex: ATTR(0, DENSITIY_SIZE, DENSITIES_OFFSET)
//  Ex: ATTR_CMP(0, HEURISTIC_SIZE, HEURISTICS_OFFSET, WORKSPACE(XMIN))
//
//=======================================/
#define ATTR(X, SIZE, OFFSET) (OFFSET + (X*SIZE))
#define ATTR_CMP(X, SIZE, OFFSET, CMP)	(OFFSET + (X*SIZE) + CMP)

// Flattens index to access RAM
#define FLATTEN_INDEX(HEATMAP, PIXEL) (HEATMAP*PIXEL)

// Map dimension components (DEPRECATED)
// #define MAPDIM_WIDTH(X) (2*X)
// #define MAPDIM_HEIGHT(X) (2*X + 1)

// Heuristic components
#define WORKSPACE(X) 		X
#define BBOX(X) 			(X + 6)
#define TOTAL_NUM_POINTS 	10

/* Auxillary defines for accessing heuristic attributes*/
#define XMIN 0
#define XMAX 1
#define YMIN 2
#define YMAX 3
#define ZMIN 4
#define ZMAX 5
/* Accessing pose components */
#define POSE_X 0
#define POSE_Y 1
#define POSE_Z 2
/* Accessing map dim components */ //(Deprecated! No longer stored in observation pyramid)
#define MAP_WIDTH 0
#define MAP_HEIGHT 1

// Heat map bounding box in each dimension
// Layer 1
#define HEATMAP_H1 9
#define HEATMAP_W1 14
// Layer 2
#define HEATMAP_H2 24
#define HEATMAP_W2 34
// Layer 3
#define HEATMAP_H3 54
#define HEATMAP_W3 74
// Layer 4
#define HEATMAP_H4 0
#define HEATMAP_W4 0
// Layer 5
#define HEATMAP_H5 0
#define HEATMAP_W5 0

// Number of bounding boxes in each heatmap
#define HEATMAP_1 (HEATMAP_H1*HEATMAP_W1)
#define HEATMAP_2 (HEATMAP_H2*HEATMAP_W2)
#define HEATMAP_3 (HEATMAP_H3*HEATMAP_W3)
#define HEATMAP_4 0
#define HEATMAP_5 0

// Memory size for sections of pyramid in words
#define MAP_DIM_SIZE  2
#define DENSITY_SIZE 	1
#define HEURISTIC_SIZE (1+4+6)

// Pyramid size by total pixels (bounding boxes in heatmap on all layers)
#define PYRAMID_PIXELS (HEATMAP_1+HEATMAP_2+HEATMAP_3+HEATMAP_4+HEATMAP_5)

// Memory offsets for pyramid
#define MAP_DIM_OFFSET 		  0      // Not stored in pyramid
#define HEURISTICS_OFFSET 	0
#define DENSITIES_OFFSET 	(HEURISTICS_OFFSET + HEURISTIC_SIZE*PYRAMID_PIXELS)

// Memory sizes for pyramid
#define HEURISTICS_MEM_SIZE (HEURISTIC_SIZE*PYRAMID_PIXELS)
#define DENSITIES_MEM_SIZE (DENSITY_SIZE*PYRAMID_PIXELS)

// Memory pyramid size in words
#define PYRAMID_MEM_SIZE  (MAP_DIM_SIZE + HEURISTIC_SIZE + DENSITY_SIZE)*PYRAMID_PIXELS

// Bounding sizes corresponding to each heatmap layer
#define BBOX_1 224*224
#define BBOX_2 112*122
#define BBOX_3 56*56
#define BBOX_4 0
#define BBOX_5 0

// Overlap amount of bounding boces
#define BBOX_OVERLAP_1 32
#define BBOX_OVERLAP_2 16
#define BBOX_OVERLAP_3 8
#define BBOX_OVERLAP_4 0
#define BBOX_OVERLAP_5 0
