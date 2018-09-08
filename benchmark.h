#pragma once
#include <vector>
#include <utils/cvl/pose.h>

/**
 * Klas p3p goes H solver for the p3p problem. Further speedup possible, comparable to kneips solver,
 * the functions and data is generated from the matlabfiles
 * The define KLAS_P3P_WITH_TIMING to test speed
 */




namespace mlib{
namespace klas{

/**
 * @brief testAll
 * tests klas sovler and compares it to kneips solver
 */
void testAll();


}// end namespace klasp3p
}// end namespace mlib

