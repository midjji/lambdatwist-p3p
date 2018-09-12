#pragma once
#include <Eigen/Dense>

#include <data.h>
//#include <mlib/utils/string_helpers.h>
//#include <mlib/utils/matlab_helpers.h>
#include "utils/random_vectors.h"


namespace cvl{



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

    Generator(){

    }

    bool degenerate_configuration(double threshold){

    }

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
            //Vector2d y((double)((float)randn(0,1)),(double)((float)randn(0,1)));
            Vector2<T> y(mlib::randu<long double>(-1,1),mlib::randu<long double>(-1,1));
          // y*=0.1;
            T z=mlib::randu<long double>(0.1,10);// 0.1 to 10 is good
            //T z=mlib::randu<long double>(0.9,1.1);// 0.1 to 1000 is good too

            Vector3<T> xri=z*y.homogeneous();
            Vector3<T> x=(P.inverse()*xri);
            x0[i]=x;
            xr[i]=xri;

           // y+=Vector2<T>(mlib::randn<long double>(0.0,0.01),mlib::randn<long double>(0,0.01));
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



};

}
