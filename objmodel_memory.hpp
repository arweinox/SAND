#pragma once

/*
 * This header files contains macros for accessing
 * and manipulating object models stored
 * within the DDR4.
 */

#include "memory.hpp"

// Number of objects
#define OBJECTS 15

// Object Model limits
#define OBJMODEL_V_MAX 1231
#define OBJMODEL_F_MAX 2000

// Count for object model vertices
#define OBJ_1V 1021
#define OBJ_2V 1002
#define OBJ_3V 1231
#define OBJ_4V 1000
#define OBJ_5V 1000
#define OBJ_6V 1002
#define OBJ_7V 1033
#define OBJ_8V 1002
#define OBJ_9V 1002
#define OBJ_10V 1002
#define OBJ_11V 1002
#define OBJ_12V 1002
#define OBJ_13V 1000
#define OBJ_14V 1010
#define OBJ_15V 1000

// Count for object model faces
#define OBJ_1F 2000
#define OBJ_2F 2000
#define OBJ_3F 1999
#define OBJ_4F 2000
#define OBJ_5F 2000
#define OBJ_6F 2000
#define OBJ_7F 1999
#define OBJ_8F 2000
#define OBJ_9F 2000
#define OBJ_10F 2000
#define OBJ_11F 2000
#define OBJ_12F 2000
#define OBJ_13F 2000
#define OBJ_14F 2000
#define OBJ_15F 2000

// Object model memory size (measured in words)
// [If the data type of faces and vertices changes such that either is no longer
// 32-bit, then you must change the bit size for each]
#define OBJ_V_WORDS 3
#define OBJ_F_WORDS 3

// Offsets for object model vertices
#define OBJ_1V_OFFSET OBJMODEL_OFFSET
#define OBJ_2V_OFFSET (OBJ_1F_OFFSET + (OBJ_F_WORDS*OBJ_1F))
#define OBJ_3V_OFFSET (OBJ_2F_OFFSET + (OBJ_F_WORDS*OBJ_2F))
#define OBJ_4V_OFFSET (OBJ_3F_OFFSET + (OBJ_F_WORDS*OBJ_3F))
#define OBJ_5V_OFFSET (OBJ_4F_OFFSET + (OBJ_F_WORDS*OBJ_4F))
#define OBJ_6V_OFFSET (OBJ_5F_OFFSET + (OBJ_F_WORDS*OBJ_5F))
#define OBJ_7V_OFFSET (OBJ_6F_OFFSET + (OBJ_F_WORDS*OBJ_6F))
#define OBJ_8V_OFFSET (OBJ_7F_OFFSET + (OBJ_F_WORDS*OBJ_7F))
#define OBJ_9V_OFFSET (OBJ_8F_OFFSET + (OBJ_F_WORDS*OBJ_8F))
#define OBJ_10V_OFFSET (OBJ_9F_OFFSET + (OBJ_F_WORDS*OBJ_9F))
#define OBJ_11V_OFFSET (OBJ_10F_OFFSET + (OBJ_F_WORDS*OBJ_10F))
#define OBJ_12V_OFFSET (OBJ_11F_OFFSET + (OBJ_F_WORDS*OBJ_11F))
#define OBJ_13V_OFFSET (OBJ_12F_OFFSET + (OBJ_F_WORDS*OBJ_12F))
#define OBJ_14V_OFFSET (OBJ_13F_OFFSET + (OBJ_F_WORDS*OBJ_13F))
#define OBJ_15V_OFFSET (OBJ_14F_OFFSET + (OBJ_F_WORDS*OBJ_14F))

// Offsets for object model faces
#define OBJ_1F_OFFSET (OBJ_1V_OFFSET + (OBJ_V_WORDS*OBJ_1V))
#define OBJ_2F_OFFSET (OBJ_2V_OFFSET + (OBJ_V_WORDS*OBJ_2V))
#define OBJ_3F_OFFSET (OBJ_3V_OFFSET + (OBJ_V_WORDS*OBJ_3V))
#define OBJ_4F_OFFSET (OBJ_4V_OFFSET + (OBJ_V_WORDS*OBJ_4V))
#define OBJ_5F_OFFSET (OBJ_5V_OFFSET + (OBJ_V_WORDS*OBJ_5V))
#define OBJ_6F_OFFSET (OBJ_6V_OFFSET + (OBJ_V_WORDS*OBJ_6V))
#define OBJ_7F_OFFSET (OBJ_7V_OFFSET + (OBJ_V_WORDS*OBJ_7V))
#define OBJ_8F_OFFSET (OBJ_8V_OFFSET + (OBJ_V_WORDS*OBJ_8V))
#define OBJ_9F_OFFSET (OBJ_9V_OFFSET + (OBJ_V_WORDS*OBJ_9V))
#define OBJ_10F_OFFSET (OBJ_10V_OFFSET + (OBJ_V_WORDS*OBJ_10V))
#define OBJ_11F_OFFSET (OBJ_11V_OFFSET + (OBJ_V_WORDS*OBJ_11V))
#define OBJ_12F_OFFSET (OBJ_12V_OFFSET + (OBJ_V_WORDS*OBJ_12V))
#define OBJ_13F_OFFSET (OBJ_13V_OFFSET + (OBJ_V_WORDS*OBJ_13V))
#define OBJ_14F_OFFSET (OBJ_14V_OFFSET + (OBJ_V_WORDS*OBJ_14V))
#define OBJ_15F_OFFSET (OBJ_15V_OFFSET + (OBJ_V_WORDS*OBJ_15V))

// Total memory size for object vertices
//#define OBJ_V_MEM_SIZE OBJ_V_WORDS*(OBJ_1V + OBJ_2V + OBJ_3V + OBJ_4V + OBJ_5V + \
//															      OBJ_6V + OBJ_7V + OBJ_8V + OBJ_9V + OBJ_10V + \
//															      OBJ_11V + OBJ_12V + OBJ_13V + OBJ_14V + OBJ_15V)

#define OBJ_V_MEM_SIZE (OBJ_V_WORDS*OBJ_1V)
// Total memory size for object faces
//#define OBJ_F_MEM_SIZE OBJ_F_WORDS*(OBJ_1F + OBJ_2F + OBJ_3F + OBJ_4F + OBJ_5F + \
//																	  OBJ_6F + OBJ_7F + OBJ_8F + OBJ_9F + OBJ_10F + \
//																	  OBJ_11F + OBJ_12F + OBJ_13F + OBJ_14F + OBJ_15F)
#define OBJ_F_MEM_SIZE (OBJ_F_WORDS*OBJ_1F)

// Total memory size for all objmodels
#define OBJMODEL_MEM_SIZE (OBJ_V_MEM_SIZE + OBJ_F_MEM_SIZE)
