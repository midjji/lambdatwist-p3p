#pragma once
/* ********************************* FILE ************************************/
/** \file    random.h
 *
 * \brief    This header contains contains convenience functions for repetable random values
 *
 * \remark
 * - c++11
 * - no dependencies
 * - self contained(just .h,.cpp)
 * - C++11 accurate generators are consistently faster and better than rand ever was!
 * - Repeatability is guaranteed if the seed is set and only these functions are used.
 *
 *
 *
 * \todo
 *
 * - convert to pure header for ease of inclusion and to ensure current not original flags are used!
 *
 *  Cmake Options: sets the flags
 * -DRANDOM_SEED_VALUE=0
 * -DRANDOM_SEED_FROM_TIME ON
 *  * option:
 *
 *
 *  * RANDOM DEFAULTS TO
 * RANDOM_SEED_VALUE 0
 * RANDOM_SEED_FROM_TIME OFF
 *
 *  Note, there is no way to have synch values for a multithreaded system.
 *
 *
 * \author   Mikael Persson
 * \date     2007-04-01
 * \note MIT licence
 *
 ******************************************************************************/
//////////////////// SELF CONTAINED ////////////////////////////




#include <random>
#include <vector>
#include <set>
#include <mutex>
#include <algorithm>
#include <assert.h>
#include <utils/cvl/pose.h>
namespace mlib{

namespace random{
// this is not thread safe, and it cant meaningfully be so either.
// no not even if you make one per thread
// no threaded program is repeatable, ever, regardless
// even in the single threaded case,
// you cannot use the generator in a static init, due to init order fiasco...

static const std::uint64_t seed{
#ifdef RANDOM_SEED_FROM_TIME
    std::chrono::system_clock::now().time_since_epoch().count()
#else
#ifndef RANDOM_SEED_VALUE
#define RANDOM_SEED_VALUE 0
#endif
RANDOM_SEED_VALUE
#endif
};
// recommended for scientific use,
// slower in some benchmarks, but faster in practice,
// especially for clang10 compared to the standard one (almost a factor of 10x faster)
static std::mt19937_64 generator(seed);
} // end namespace random

double randu(double low=0, double high=1);
int randui(int low=0, int high=1);
double randn(double mean=0, double sigma=1);


} // end namespace mlib
namespace cvl {


PoseD random_pose();
Matrix3d random_rotation_matrix();
}
