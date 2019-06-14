#pragma once //#ifndef CONFIG

/* Header contains all necessary settings for algorithm.
 * For detailed, settings regarding observation pyramid
 * or object models edit the appropriate header file
 */

/* SIMULATION****************************************/
//#define SIMULATION // Disable for none HLS code, enable for simulation feedback

/* CAMERA VIEW**************************************/
// Indicates if 3D model is rendered with same camera view of sensor
#define SAME_CAMERAVIEW
/* Algorithm options*******************************/
#define DIFFERENT_BBOXES
//#define USE_NORMAL
//#define ENABLE_VISUALIZATION
#define NUM_OBJECTS 1
/* PYRAMID SETTINGS*********************************/
// Pyramid size by heatmap layers
#define PYRAMID_HEATMAPS     3
#define MAX_PYRAMID_HEATMAPS 5
/* IMAGE DIMENSIONS*********************************/
#define WIDTH 480
#define HEIGHT 640
#define IMG_SIZE (WIDTH*HEIGHT)
/* CAMERA SETTINGS**********************************/
// Camera properties
#define NEAR 0.1
#define FAR 1000
// Camera intrinsic propteries
#define FX 542.556
#define FY 542.566
#define CX 332.065
#define CY 222.214

// Camera view Matrix
//  first row
#define UX 0.00793759
#define UY -0.999892
#define UZ 0.0123589
#define UW 0.00663925
//  second row
#define VX 0.787717
#define VY 0.0138718
#define VZ 0.615881
#define VW -0.954239
//  third row
#define NX -0.615986
#define NY 0.00485456
#define NZ 0.787742
#define NW -0.976851

/* IterativeLikelihoodReweighting configuration parameters*/
// Something to do with managing pixels
#define HORIZONTAL_COUNT 25
#define VERTICAL_COUNT 25
// Number of samples
#define NUM_SAMPLES 10
// Number of iterations
#define NUM_ITERATIONS 1
#define ENABLE_SIGMOID
#define DECAY 0.0001
#define TRANSITION_VARIANCE 0.08
#define ROTATION_VARIANCE 0.4
#define CONVERGENCE_THRESHOLD 0.7
#define NORM_RANGE_FOR_SIGMOID 3
// Iteration parameters
#define RESTART_PER_ITERATION 50
#define NUM_RESTART_SAMPLE 625
#define MAX_RESTART_TIMES 4
// Coefficient values
#define INLIERS_IN_BBOX 0.3
#define INLIERS_IN_RENDER 0.7
// Table height (represents minimum z for workspace)
#define TABLE_HEIGHT 0.75

//#endif /* end of include guard: CONFIG */
