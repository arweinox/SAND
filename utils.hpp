#pragma once
//
// @file utils/utils.h
#include "sandtypes.hpp"

#ifndef __SYNTHESIS__
#include <stdio.h>
#include "config.hpp"
#include "memory.hpp"
#include "sample.hpp"
#endif


// @brief Utils
// @author Kevin Anderson


/* Definition from the GLM library */
#define PI 3.14159265358979323846264338327950288

// May need to change from pass-by-reference
void RoundAngle(float& radian);

/* Workaround for pointer reinterperation (int to float)
 * from dram bus
 *
 * dram - bus to dram from top function
 * index - location in dram for data
 * return data reinterpreted as float
 */
float ReinterpretFloat(in_32* dram, int index);

float BoundValue(float val, float lbound, float ubound);

/*
 * Helper function to print out relevant information for
 * each sample
 */
void PrintSamples(in_32* dram);
