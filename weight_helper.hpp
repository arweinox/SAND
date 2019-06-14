#pragma once
//
// @file utils/utils.h
// @brief Utils
// @author Zhiqiang Sui
// University of Michigan, 2017
//

#include <cmath>

inline float Sigmoid(const float& val)
{
    return 1 / (1 + exp(-val));
}
