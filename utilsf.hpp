// Jan 25, 2019: Purpose is to contain helper functions done in cuda for inlier computation

#pragma #ifndef UTILS
#define UTILS

#define ROWS height*width
#define COLS horizontal_count*vertical_count

static float *objective_buffer;
// weight for each sample
static float weights[COLS];

// For computing the average particles
static float *allparticles;
static float average[ROWS];
static float *avgdivisor;

// For computing number of points rendered in the bounding box
static float num_rendered_points[COLS];
static float num_rendered_points_in_bbox[COLS];
static float *rendered_points;
static float *rendered_points_in_bbox;


#endif /* end of include guard: UTILS */
