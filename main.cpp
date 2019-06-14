#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include "sample.hpp"
//#include "config.hpp"
//#include "sandtypes.hpp"
#include "utils.hpp"
#include "sandtop.hpp"
#include "bluecup_densities.hpp"
#include "bluecup_heuristics.hpp"
#include "bluecup_objmodelf.hpp"
#include "bluecup_objmodelv.hpp"
#include "scene_exp01.hpp"
#include "objmodel_memory.hpp"

using namespace std;


// March 13, 2019: This main function does no testing. Only issues dummy call to
// reveal build errors

int main() {

	/******************** INITIALIZE VARIABLES FOR FUNCTION TO BE TESTED */
	int allocate = PYRAMID_MEM_SIZE +
					OBJMODEL_MEM_SIZE +
					WIDTH*HEIGHT +
					NUM_SAMPLES*RESULTS_MEM_SIZE + 5000;

	in_32 * dram = (in_32 *) malloc(allocate*sizeof(int));

	// Set pyramid heuristics
	for(int y = 0; y <  HEURISTICS_MEM_SIZE; y++) {
		dram[y + HEURISTICS_OFFSET] = bluecup_heuristics[y];
	}
	// Set pyramid densities
	for(int y = 0; y <  DENSITIES_MEM_SIZE; y++) {
		dram[y + DENSITIES_OFFSET] = bluecup_densities[y];
	}
	// Set observation
	for(int x = 0; x < IMG_SIZE; x++) {
		dram[OBSERVATION_OFFSET + x] = reinterpret_cast<int &>(scene_exp01[x]);
	}
	// Set model vertices
	for(int x = 0; x < OBJ_V_MEM_SIZE; x++) {
		dram[OBJMODEL_OFFSET  +  x] = reinterpret_cast<int &>(bluecup_objmodelv[x]);
	}
	// Set model face
	for(int x = 0; x < OBJ_F_MEM_SIZE; x++) {
		dram[OBJMODEL_OFFSET + OBJ_V_MEM_SIZE + x] = bluecup_objmodelf[x];
	}



	/*********************************************************************/

	// Initialize return variable
	int ret = 0;
	ap_uint<1> done[] = {0};
	ap_uint<STATUS_STATES> status[] = {0};

	// Make function call
//	sandtop(dram, done, status);

	// Write weights to file
//	PrintSamples(dram);

	// Verify results
//	ret = system("diff --brief -w output.dat correct_output.dat");
	if(ret) {
		printf("\nTest failed\n\n");
	} else {
		printf("\nTest passed\n\n");
	}

	// Return 1 if test failure, 0 otherwise
	return ret;
}
