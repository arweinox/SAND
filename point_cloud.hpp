#pragam #ifndef POINTCLOUD
#define POINTCLOUD

// ************DEPRACATED!*********************/

#include "config.hpp"

typedef float[weight][height] point_cloud;

struct point_cloud {
  // May not need internal sizes; sizes should be declared in config.hpp
  int width;
  int height;
  point_cloud pc;
}

#endif /* end of include guard: POINTCLOUD */
