#pragma once


//'terms' is defined as std::map<string, float> terms_; in the orignal C++ code
//For this implementation, one sample has a single object,

//4 32-bit words
struct Weight
{
	// Coef
	float terms_[2];

	float final_weight_;
};
