#pragma once
#include <utils/mlibtime.h>
#include <data.h>

class P3PResult{
public:
    P3PResult(std::string name, int capacity):name(name){errors.reserve(capacity);}

    std::string name;
    // number of times that not a single valid solution was output!
    int no_solution=0;

    // number of times the ground truth was found
    int ground_truth_in_set=0;

    // correct solutions according to solver
    int valid=0;

    // incorrect solutions output by solver
    int incorrect_valid=0;

    // number of unique solutions <= valid
    int solutions=0;

    // number of correct duplicates
    int duplicates=0;
    // vector with the smallest error, or 1 for nan etc...
    std::vector<float128> errors;
};






