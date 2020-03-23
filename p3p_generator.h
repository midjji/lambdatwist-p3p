#pragma once
//#include <Eigen/Dense>

#include <data.h>

#include "utils/random_vectors.h"


namespace cvl{


/*
double angle(Vector3d a, Vector3d b){
    a.normalize();
    b.normalize();
    return std::acos(a.dot(b));
}
double angle(Vector2d a, Vector2d b){
    a.normalize();
    b.normalize();
    return std::acos(a.dot(b));
}
*/


/**
 * @brief The Generator class
 *
 * generates samples for the p3p problem,
 * uniform rotation sampling by the projected gaussian trick,
 * gaussian translations since the scale of the translation can always be eliminated from the problem anyways.
 *
 * picks image coordinates randomly from a pinhole normalized camera spanning -1,1. always possible to set.
 */
class Generator{
public:

    Generator(){}


    template<class T>
    Data<T> next(){

        cvl::Matrix<T,3,3> R=mlib::getRandomRotation<double>();
        cvl::Vector3<T> t(mlib::randn<long double>(0,1),mlib::randn<long double>(0,1),mlib::randn<long double>(0,1));
        t.normalize();
        Pose<T> P(R,t);
        Vector3<Vector3<T>>  x0,xr;
        Vector3<Vector2<T>>  ys;



        // good generator, pick a uniform position between -1,1 and multiply by a distance...
        for(int i=0;i<3;++i){
            Vector2<T> y(mlib::randu<long double>(-1,1),mlib::randu<long double>(-1,1));
            T z=mlib::randu<long double>(0.1,100);

            Vector3<T> xri=z*y.homogeneous();
            Vector3<T> x=(P.inverse()*xri);
            x0[i]=x;
            xr[i]=xri;
            ys[i]=y;
        }


        // collinearity or identity are degenerate configurations and should be avoided.
        bool colinear=false;
        {
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    if(i!=j){
                        if(x0(i)==x0(j) || ys[i]==ys[j])
                            colinear=true;
                    }
        }
        if(colinear) return next<T>();
        return Data<T>(P,x0,xr);
    }


    // there are a number of special cases
    // note that many are symmetric, but all should be included because the solvers might not be.
    std::vector<Data<T>> special_cases(){
        // these are the hard special cases, ie exact

        // special case 1:3 are
        // y_i'*y_j = 0 for 12, 13, 23
        // special case 4:6 are
        // |x_i - x_j| = 0
        // special case
        //







    }




};

}
