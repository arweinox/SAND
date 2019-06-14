 #pragma once

#include "../source/config.hpp"
#include "../source/object.hpp"
#include "../source/sandtypes.hpp"
#include "../source/utils.hpp"
#include "global.hpp"

#ifndef __SYNTHESIS__
  #include <cmath>
  #include <string>
  #include <cstdio>
  #include <limits>
  #include <iostream>
  #include <fstream>
#include <stdio.h>
#endif

#define DIM 4
// #define LEFT -CX
// #define RIGHT (WIDTH + CX)
// #define TOP (HEIGHT + CY)
// #define BOTTOM -CY
// Edge function expressed as cross product. More efficient than expanding form
// where C term has 2 multiplies which totaling to 4 multiplies and 3 additions
#define EDGE(X0,Y0,X1,Y1,X2,Y2) ((X2-X0)*(Y1-Y0))-((Y2-Y0)*(X1-X0))  // X0  - the common point, X1 - another point, X2 - is the point you want to check if it is inside the triangle
#define SCALENDC(X, ULIMIT) (ULIMIT*(X + 1)/2)
#define SCALENDCDEPTH(D) ((-FAR+NEAR)*(D + 1)/2 - NEAR)   // [-1,1] - > [-NEAR, -FAR]

// Represent triangle for rendering
typedef sandy::int3 triangle;

namespace renderer {
template <class T>
class Renderer {
public:
  Renderer() {
    for(int i = 0; i < IMG_SIZE; i++) {
      z_buffer_[i] = 0;
    }
  }

  ~Renderer() {}

  /*
   * Initialize the matrices with elements
   *
   * sample_pose - pose for current sample
   * dim - x,y,z scaling factors for current sample
   * near - variable needed
   * far - variable needed
   */
  void InitializeMatricies(int sample_index) {
#ifndef __SYNTHESIS__
#ifdef VERBOSE
			std::cout << "Initialize Matrices (Renderer)" << std::endl;
#endif
#endif

    // Extract euler angles
    float ex = sampling::sample_list_[sample_index].object_list_.pose_.euler_[0];
    float ey = sampling::sample_list_[sample_index].object_list_.pose_.euler_[1];
    float ez = sampling::sample_list_[sample_index].object_list_.pose_.euler_[2];

    // Extract cartesian coordinates
    float X = sampling::sample_list_[sample_index].object_list_.pose_.pos_[0];
    float Y = sampling::sample_list_[sample_index].object_list_.pose_.pos_[1];
    float Z = sampling::sample_list_[sample_index].object_list_.pose_.pos_[2];

    // Extract scale factor
    float sx = sampling::sample_list_[sample_index].object_list_.dim_[0];
    float sy = sampling::sample_list_[sample_index].object_list_.dim_[1];
    float sz = sampling::sample_list_[sample_index].object_list_.dim_[2];

    // Initialize scaling matrix
    S[0][0] = sx;
    S[1][1] = sy;
    S[2][2] = sz;
    S[3][3] = 1;
    for(int i = 0; i < DIM; i++) {
      for(int j = 0; j < DIM; j++) {
        if(i != j){ S[i][j] = 0; }
      }
    }

    // Initialize translation matrix
    Tr[0][3] = X;
    Tr[1][3] = Y;
    Tr[2][3] = Z;
    Tr[3][3] = 1;
    for(int i = 0; i < DIM; i++) {
      for(int j = 0; j < DIM-1; j++) {
        i==j ? Tr[i][j] = 1 : Tr[i][j] = 0;
      }
    }

    // Initialize x-axis rotation matrix
    Rx[0][0] = 1;
    Rx[0][1] = 0;
    Rx[0][2] = 0;
    Rx[0][3] = 0;

    Rx[1][0] = 0;
    Rx[1][1] = cos(ex);
    Rx[1][2] = -sin(ex);
    Rx[1][3] = 0;

    Rx[2][0] = 0;
    Rx[2][1] = sin(ex);
    Rx[2][2] = cos(ex);
    Rx[2][3] = 0;

    Rx[3][0] = 0;
    Rx[3][1] = 0;
    Rx[3][2] = 0;
    Rx[3][3] = 1;

    // Initialize y-axis rotation matrix
    Ry[0][0] = cos(ey);
    Ry[0][1] = 0;
    Ry[0][2] = sin(ey);

    Ry[0][3] = 0;

    Ry[1][0] = 0;
    Ry[1][1] = 1;
    Ry[1][2] = 0;
    Ry[1][3] = 0;

    Ry[2][0] = -sin(ey);
    Ry[2][1] = 0;
    Ry[2][2] = cos(ey);
    Ry[2][3] = 0;

    Ry[3][0] = 0;
    Ry[3][1] = 0;
    Ry[3][2] = 0;
    Ry[3][3] = 1;

    // Initialize z-axis rotation matrix
    Rz[0][0] = cos(ez);
    Rz[0][1] = -sin(ez);
    Rz[0][2] = 0;
    Rz[0][3] = 0;

    Rz[1][0] = sin(ez);
    Rz[1][1] = cos(ez);
    Rz[1][2] = 0;
    Rz[1][3] = 0;

    Rz[2][0] = 0;
    Rz[2][1] = 0;
    Rz[2][2] = 1;
    Rz[2][3] = 0;

    Rz[3][0] = 0;
    Rz[3][1] = 0;
    Rz[3][2] = 0;
    Rz[3][3] = 1;

    // Initialize camera view matrix
    C[0][0] = UX;
    C[0][1] = UY;
    C[0][2] = UZ;
    C[0][3] = UW;

    C[1][0] = VX;
    C[1][1] = VY;
    C[1][2] = VZ;
    C[1][3] = VW;

    C[2][0] = NX;
    C[2][1] = NY;
    C[2][2] = NZ;
    C[2][3] = NW;

    C[3][0] = 0;
    C[3][1] = 0;
    C[3][2] = 0;
    C[3][3] = 1;

    // Initialize projection matrix
    P[0][0] = 2 * FX / WIDTH;
    P[0][1] = 0;
    P[0][2] = 1 - (2 * CX / WIDTH);
    P[0][3] = 0;

    P[1][0] = 0;
    P[1][1] = 2 * FY / HEIGHT;
    P[1][2] = 1 - (2 * CY / HEIGHT);
    P[1][3] = 0;

    P[2][0] = 0;
    P[2][1] = 0;
    P[2][2] = -(FAR + NEAR)/(FAR - NEAR);
    P[2][3] = -(2 * FAR * NEAR) / (FAR - NEAR);

    P[3][0] = 0;
    P[3][1] = 0;
    P[3][2] = -1;
    P[3][3] = 0;
  }

  /*
   * Performs matrix multiplication for square matricies. Based on Xilinx
   * Application Notes: A Zynq Accelerator for Floating Point
   * Matrix Multiplication Designed with
   * Vivado HLS
   *
   * A - 4x4 matrix, left-hand matrix for multiplication
   * B - 4x4 matrix, input coordinates
   * C - 4x4 matrix, tranformed coordinates
   */
  void MatrixMult(T A[DIM][DIM], T B[DIM][DIM], T C[DIM][DIM]) {
    // matrix multiplication of a A*B matrix
    for (int ia = 0; ia < DIM; ++ia) {
      for (int ib = 0; ib < DIM; ++ib) {
        T sum = 0;
        for (int id = 0; id < DIM; ++id) {
          sum += A[ia][id] * B[id][ib];
        }
        C[ia][ib] = sum;
      }
    }
  }

  /*
   * Performs matrix multiplication for square matrices. Based on Xilinx
   * Application Notes: A Zynq Accelerator for Floating Point
   * Matrix Multiplication Designed with
   * Vivado HLS
   *
   * A - 4x4 matrix, left-hand matrix for multiplication
   * B - 4x1 matrix, input coordinates
   * C - 4x1 matrix, transformed coordinates
   */
  void TransformMatrixMult(sandy::float3 point, float* point_t) {
    // matrix multiplication of a A*B matrix
    point_t[0] = FT[0][0]*point.x + FT[0][1]*point.y + FT[0][2]*point.z + FT[0][3];
    point_t[1] = FT[1][0]*point.x + FT[1][1]*point.y + FT[1][2]*point.z + FT[1][3];
    point_t[2] = FT[2][0]*point.x + FT[2][1]*point.y + FT[2][2]*point.z + FT[2][3];
    point_t[3] = FT[3][0]*point.x + FT[3][1]*point.y + FT[3][2]*point.z + FT[3][3];
  }

  /*
   * This function performs rigid-body transformations and transforms points
   * into camera viewpoint
   *
   * object_pc_ - bus to RAM holding 3D object model
   * num_points - the number of vertices in model
   */
  void Transform(int num_points) {    // Needs to take point cloud BRAM as argument
#ifndef __SYNTHESIS__
	#ifdef VERBOSE
			std::cout << "Transform (Renderer)" << std::endl;
	#endif
#endif

    T temp[DIM][DIM];
    T temp2[DIM][DIM];
//    T point[DIM][1];
    T point_t[DIM];    // transformed point
    // Create rigid-body transformation Matrix (Rx and temp temporarily holds results for next multiplication)
    MatrixMult(Ry, Rx, temp);
    MatrixMult(Rz, temp, temp2);
    MatrixMult(Tr, temp2, temp);
    MatrixMult(S, temp, temp2);
    MatrixMult(C, temp2, CSTR);
    MatrixMult(P, CSTR, FT);      // Matrix embodying rigid-body-transform, camera view and projection matrix

    // Transform entire object. Iterate over entire point cloud of object.
    for(int point = 0; point < num_points; point++) {
      TransformMatrixMult(object_modelv_[point], point_t);
      // Normalized by w
      object_modelvt_[point].x = point_t[0]/point_t[3];
      object_modelvt_[point].y = point_t[1]/point_t[3];
      object_modelvt_[point].z = point_t[2]/point_t[3];   // Do I need to divide to to get into NDC if I wnat the depth
    }
  }

  /*
   * This function performs rasterisation on 3D object model constructed from triangles
   * of known vertices. The entire bounding box for the object is search, however,
   * for efficiency, multiplications involved in the edge function are substituted
   * with simple comparisons and additions. Edge functions are used to tell whether
   * a give pixel of the resulting depth map lies within triangle. These functions
   * compute the magnitude of the cross product between vectors form be an edge
   * and the propose point with a common vertex.
   *
   * faces - RAM storing the triangles which hold the indices of vertices in ascending order
   * vertices - RAM storing vertices sorted by x-coordinate in ascending order
   * z-buffer - RAM representing depth map
   * face_num - number of faces in model
   */
   /*
    * Possible optimizations:
    *   change for-loops to while loops
    *   reduce bits of i and j
    *   reduce bits for all precomputes value with fixed point
    */
  void Rasterize(int face_num) {
#ifndef __SYNTHESIS__
	#ifdef VERBOSE
				std::cout << "Rasterize (Renderer)" << std::endl;
	#endif
#endif

//    #ifndef __SYNTHESIS__
//	  FILE * file;
//	  file = fopen("rasterized.txt", "w+");
//      if(!file) {
//        std::cout << "failed to open rasterized" << std::endl;
//      }
//    #endif

    int i, j;                     // pixels indices, row and columns respectively
    int offset;                   // column offset for access z-buffer for raster order
    float area;										// Area of face (triangle)
    float A12, A23, A31;          // x-coordinate  coefficient in the expanding cross product
    float B12, B23, B31;          // y-coordinate  coefficient in the expanding cross product
    float alpha, beta, gamma;     // Barycenteric coordinates
    renderer::BoundingBox2D box;  // bounding box to iterate over for face

    int total_count = 0;
    for(int x = 0; x < face_num; x++) { // iterate through each face
       bool raster_flag = false;          // Flag for searching in raster order
       // Access vertices for current face
       triangle face = object_modelf_[x];
#ifndef __SYNTHESIS__
//       std::cout << "x: " << x << std::endl;
//       printf("		face a: %d, face b: %d, face c: %d\n", face.a, face.b, face.c);
       if(face.a < 0 || face.a > 1020 || face.b < 0 || face.b > 1020 || face.c < 0 || face.c > 1020) {
    	   std::cout << " Problem" << std::endl;
       }
#endif
       sandy::float3 v1 = object_modelvt_[face.a];
       sandy::float3 v2 = object_modelvt_[face.b];
       sandy::float3 v3 = object_modelvt_[face.c];
//#ifndef __SYNTHESIS__
//       std::cout << "	face loaded" << std::endl;
////     			printf("face a: %d, face b: %d, face c: %d\n", face.a-1, face.b-1, face.c-1);
//#endif

       // Scale NDC to image size to convert to pixels
       v1.x = SCALENDC(v1.x, WIDTH);   // Might be WIDTH
       v2.x = SCALENDC(v2.x, WIDTH);   // Might be WIDTH
       v3.x = SCALENDC(v3.x, WIDTH);   // Might be WIDTH

       v1.y = SCALENDC(v1.y, HEIGHT);    // Might be HEIGHT
       v2.y = SCALENDC(v2.y, HEIGHT);    // Might be HEIGHT
       v3.y = SCALENDC(v3.y, HEIGHT);    // Might be HEIGHT

       // v1.z = SCALENDCDEPTH(v1.z);    // Might be HEIGHT
       // v2.z = SCALENDCDEPTH(v2.z);    // Might be HEIGHT
       // v3.z = SCALENDCDEPTH(v3.z);    // Might be HEIGHT
       //
       // Coefficients for edge equations
       A12 = v2.y - v1.y;
       A23 = v3.y - v2.y;
       A31 = v1.y - v3.y;
       B12 = v1.x - v2.x;
       B23 = v2.x - v3.x;
       B31 = v3.x - v1.x;
       // Area of triangle
       area = EDGE(v1.x, v1.y, v2.x, v2.y, v3.x, v3.y);
       // Will probably implement these are regs and not wires
       box.x_min_ = v1.x < v2.x ? (v1.x < v3.x ? v1.x : v3.x) : (v2.x < v3.x ? v2.x : v3.x);
       box.x_max_ = v1.x > v2.x ? (v1.x > v3.x ? v1.x : v3.x) : (v2.x > v3.x ? v2.x : v3.x);
       box.y_min_ = v1.y < v2.y ? (v1.y < v3.y ? v1.y : v3.y) : (v2.y < v3.y ? v2.y : v3.y);
       box.y_max_ = v1.y > v2.y ? (v1.y > v3.y ? v1.y : v3.y) : (v2.y > v3.y ? v2.y : v3.y);
       // Line equations for vertices (start upper-left corner of bounding box)
       float E12 = EDGE(v1.x, v1.y, v2.x, v2.y, box.x_min_, box.y_min_);
       float E23 = EDGE(v2.x, v2.y, v3.x, v3.y, box.x_min_, box.y_min_);
       float E31 = EDGE(v3.x, v3.y, v1.x, v1.y, box.x_min_, box.y_min_);
       // Limit bounding box to image size
       box.x_min_ = BoundValue(box.x_min_, 0, WIDTH);
       box.x_max_ = BoundValue(box.x_max_, 0, WIDTH);
       box.y_min_ = BoundValue(box.y_min_, 0, HEIGHT);
       box.y_max_ = BoundValue(box.y_max_, 0, HEIGHT);

       /////////////////////////////////////////////////
       if(v2.y < 0) {
    	   v2.y = -v2.y;
       }
       if(v1.y < 0) {
    	   v1.y = -v1.y;
		}
       if(v3.y < 0) {
     	   v3.y = -v3.y;
 		}
       /////////////////////////////////////////////////

       if(v2.y < v1.y) {    // Negate cross products (reverse counter-clockwise winding; triangle is upside down)
         // Reorient cross-product
         E12 = -E12;
         E23 = -E23;
         E31 = -E31;
         area = -area;
         // Negate sign of increments to agree with
         // orientation of cross product
         A12 = -A12;
         A23 = -A23;
         A31 = -A31;
         B12 = -B12;
         B23 = -B23;
         B31 = -B31;
       }

       // Floor bounding box dimensions for pixels
       box.x_min_ = floorf(box.x_min_);
       box.x_max_ = ceilf(box.x_max_);
       box.y_min_ = floorf(box.y_min_);
       box.y_max_ = ceilf(box.y_max_);

       /******************************************/
//       #ifndef __SYNTHESIS__
//			 	  int count = 0;
//          fprintf(file, "face: %d\n", x);
//          fprintf(file, "  Face: [%d %d %d]\n",  face.a, face.b, face.c);
//          fprintf(file, "  Vertices: [%f %f; %f %f; %f %f]\n",  reinterpret_cast<float &>(v1.x), v1.y, v2.x, v2.y, v3.x, v3.y);
//          fprintf(file, "  Bounding box: [%f %f; %f %f]\n",  box.x_min_, box.x_max_, box.y_min_, box.y_max_);
//          fprintf(file, "  E12: %f\n", E12);
//          fprintf(file, "  E23: %f\n", E23);
//          fprintf(file, "  E31: %f\n", E31);
//          fprintf(file, "  area: %f\n", area);
//          if(v2.y < v1.y)
//            fprintf(file, "  reorienting cross product...\n");
//       #endif

       /*****************************************/

       // i,j are row and column respectively
       for(j = box.y_min_; j <= box.y_max_; j++) {
         for(i = box.x_min_; i <= box.x_max_; i++) {
        	// printf("E12: %f 	E23: %f 	E31: %f\n", E12, E23, E31);
           offset = raster_flag == false ? i : (box.x_min_ + box.x_max_ - i);
           if((E12 >= 0.0) && (E23 >= 0.0) && (E31 >= 0.0)) {
             // Fill Z-Buffer for i,j
             alpha = E23/area;
             beta  = E31/area;
             gamma = E12/area;
             float z = alpha*v1.z + beta*v2.z + gamma*v3.z;   //
//				 #ifndef __SYNTHESIS__
//					 if(count==1) {
//					   fprintf(file, "    alpha: %f\n    beta: %f\n    gamma: %f\n    sum: %f\n", alpha, beta, gamma,  alpha + beta + gamma);
//					   fprintf(file, "    depth: %f\n", z);
//					 }
//					 count++;
//				 #endif
             // fprintf(file, "  depth(%d, %d): %f \n", j, i, z);
             z_buffer_[j*WIDTH + offset] = z < z_buffer_[j*WIDTH + offset] ? z : z_buffer_[j*WIDTH + offset];

           }
           if(i != box.x_max_) {  // Do not increment on last element in scanline
             // Increment edge equations with x-coordinate coefficient
             if(raster_flag == false){
               E12 += A12;
               E23 += A23;
               E31 += A31;
             } else {
               E12 -= A12;
               E23 -= A23;
               E31 -= A31;
             }
           }
//           #ifndef __SYNTHESIS__
//            if((E12 >= 0.0) && (E23 >= 0.0) && (E31 >= 0.0))
//              fprintf(file, "inside\n");
//            else
//              fprintf(file, "outside\n");
//           #endif
         }
         // Increment edge equations with y-coordinate coefficient
         E12 += B12;
         E23 += B23;
         E31 += B31;
         // Flip raster order
         raster_flag = !raster_flag;
       }
//       #ifndef __SYNTHESIS__
//          fprintf(file, "     points rendered: %d\n", count);
//          total_count += count;
//        #endif
    }

//	  #ifndef __SYNTHESIS__
//      printf("total rendered points: %d\n", total_count);
//      fclose(file);
//    #endif
  }

  /*
   * Returns the depth at pixel (x,y)
   *
   * x - x-coordinate of pixel
   * y - y-coordinate of pixel
   */
  float fetchdepth(in_32 x, in_32 y) {
	  return z_buffer_[x*WIDTH + y];
  }

  /*
   * Resets the entire z buffer to Kinect limit
   */
  void ResetZBuffer() {
	  for(int index = 0; index < IMG_SIZE; index++) {
#pragma HLS unroll factor=2
#pragma HLS pipeline II=1
		  z_buffer_[index] = FAR;
	  }
  }

  /*
   * Resets the entire z buffer to Kinect limit
   */
  void ResetElement(int y, int x) {
	z_buffer_[x*WIDTH + y] = FAR;
  }

  /*
   * FOR TESTING ONLY; Return pointer to depth map
   */
   float * getbuffer(){
     return z_buffer_;
   }

private:
  T Rx[DIM][DIM];   // x-axis rotation matrix
  T Ry[DIM][DIM];   // y-axis rotation matrix
  T Rz[DIM][DIM];   // z-axis rotation matrix
  T S[DIM][DIM];    // scaling matrix
  T Tr[DIM][DIM];   // translation matrix
  T C[DIM][DIM];    // camera view matrix
  T P[DIM][DIM];    // projection matrix
  T CSTR[DIM][DIM]; // rigid-body and camera viewpoint matrix
  T FT[DIM][DIM];   // Final transformation matrix; combing all matrices above

   // vertices of transformed 3D model of the current object
   sandy::float3 object_modelvt_[OBJMODEL_V_MAX];

  // Z-buffer of rendered object
  float z_buffer_[IMG_SIZE];
};
}
