CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(lambdatwist CXX) # defines ${CMAKE_PROJECT_NAME}(top level) PROJECT_NAME(most recent)

macro(print_list name list)
message("\n Listing: ${name}\n")
foreach(item IN LISTS ${list})
message("     ${name}: ${item}")
endforeach()
message("\n Listing: ${name} - done \n")
endmacro()



message("To exactly replicate the results from the paper, you must use ubuntu 14.04 and its default compiler")






if(NOT CMAKE_BUILD_TYPE)
    MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Debug Release"  FORCE)
else()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        MESSAGE("\n${line}")
        MESSAGE("\n-- Build type: Debug. Performance will be terrible!")
        MESSAGE("-- Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build.")
        MESSAGE("\n${line}")
    endif()
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-fno-rtti -fno-omit-frame-pointer -pg  -rdynamic ") # dynamic is for the improved asserts
set(CMAKE_CXX_FLAGS_RELEASE "-fno-rtti  -O3 -march=native -mtune=native -DNDEBUG")





option(RANDOM_SEED_FROM_TIME "generate a seed from time, default" OFF)
if(RANDOM_SEED_FROM_TIME)
    add_definitions(-DRANDOM_SEED_FROM_TIME)
else()
    add_definitions(-DRANDOM_SEED_VALUE=${RANDOM_SEED})
endif()




add_definitions("-std=c++11")
include_directories(".")
include_directories("utils")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
FIND_PACKAGE(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})
message("\n -- Found Eigen Version: ${Eigen_VERSION}")
message("\n -- Eigen_INCLUDE_DIRS: ${EIGEN_INCLUDE_DIRS}")


#including headers here make them show up in older ides

set(lambdatwist-srcs "")
list(APPEND lambdatwist-srcs lambdatwist/p3p_timers.h)
list(APPEND lambdatwist-srcs lambdatwist/lambdatwist.p3p.h)
list(APPEND lambdatwist-srcs lambdatwist/solve_cubic.h)
list(APPEND lambdatwist-srcs lambdatwist/solve_eig0.h)
list(APPEND lambdatwist-srcs lambdatwist/refine_lambda.h)

set(kneip-srcs "")
list(APPEND kneip-srcs kneip/kneip.h)
list(APPEND kneip-srcs kneip/kneip.cpp)
list(APPEND kneip-srcs kneip/kneip_utils.h)

set(ke-srcs "")
list(APPEND ke-srcs ke/ke.h)
list(APPEND ke-srcs ke/ke.cpp)
list(APPEND ke-srcs ke/ke_utils.h)

set(util-srcs "")

list(APPEND util-srcs utils/mlibtime.cpp)
list(APPEND util-srcs utils/mlibtime.h )
list(APPEND util-srcs  utils/matlab_helpers.h)
list(APPEND util-srcs utils/random.h)
list(APPEND util-srcs  utils/random_vectors.h)

list(APPEND util-srcs  utils/string_helpers.h)
list(APPEND util-srcs  utils/vector_utils.h)
list(APPEND util-srcs  utils/cvl/matrix.h)
list(APPEND util-srcs  utils/cvl/pose.h)
list(APPEND util-srcs utils/cvl/rotation_helpers.h)


print_list("utils" util-srcs)
#file(GLOB_RECURSE utils2 )

list(APPEND util-srcs benchmark.h)
list(APPEND util-srcs benchmark.cpp)
list(APPEND util-srcs data.h)
list(APPEND util-srcs solver.h)
list(APPEND util-srcs mnumerics.h)
list(APPEND util-srcs p3p_generator.h)

# headers are included for ide benefit...
add_executable(main main.cpp ${util-srcs} ${ke-srcs} ${kneip-srcs} ${lambdatwist-srcs})

message("================================================================================")
message("The p3p solvers are under their respective folders,")
message("Timing varies between compilers and processors, rankings are typically consistent however" )
message("By default a reduced set of tests is performed for the solvers, \n this benefits ke and kneip in terms of speed, but also means lots of bad solutions are passed on")
message("use the flag -DFULL_SOLUTION_TESTS to activate the full set of tests when performing timing comparisions.")
message("\n the compiler is alloed to use sse intrinsics by default, this improves speed at the cost of numerical accuracy.\n some compilers succeed, others do not\n
This has a significant impact on the accuracy/speed of kes method")
