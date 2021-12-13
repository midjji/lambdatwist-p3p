#pragma once
//#include <Eigen/Dense>

#include <data.h>

#include "utils/random.h"


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

    template<class T>
    Data<T> next(){

        cvl::Matrix<T,3,3> R=cvl::random_rotation_matrix();
        cvl::Vector3<T> t(mlib::randn(0,1),mlib::randn(0,1),mlib::randn(0,1));
        t.normalize();
        Pose<T> P(R,t);
        Vector3<Vector3<T>>  x0,xr;
        Vector3<Vector2<T>>  ys;



        // good generator, pick a uniform position between -1,1 and multiply by a distance...
        for(int i=0;i<3;++i){
            Vector2<T> y(mlib::randu(-1,1),mlib::randu(-1,1));
            T z=mlib::randu(0.1,100);

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

    bool first_special_case0=true;
    template<class T> Data<T> special_case0(){
        // 3d point differences are orthogonal
        // there is a faster lambdatwist special case solver for this, but it works outright
        Vector3d x0(0,0,0);
        Vector3d x1(1,0,0);
        Vector3d x2(0,1,0);
        if(!first_special_case0){
            PoseD P(random_rotation_matrix());
            x0=P*x0;
            x1=P*x1;
            x2=P*x2;
        }
        first_special_case0=false;

        Vector3d x0c,x1c,x2c;
        PoseD Pcam;
        //int tries=0;
        while(true){
          //  tries++;

            Pcam=random_pose();
            x0c=Pcam*x0;
            x1c=Pcam*x1;
            x2c=Pcam*x2;
            if(x0c[2]<1e-3 || x1c[2]<1e-3|| x2c[2]<1e-3) continue;
            Vector2d y0=x0c.dehom();
            Vector2d y1=x1c.dehom();
            Vector2d y2=x2c.dehom();
            if(!y0.is_in(Vector2d(-1,-1),Vector2d(1,1))) continue;
            if(!y1.is_in(Vector2d(-1,-1),Vector2d(1,1))) continue;
            if(!y2.is_in(Vector2d(-1,-1),Vector2d(1,1))) continue;
            break;
        }
        //std::cout<<"tries: "<<tries<<std::endl;// takes about 100x each

        return Data<T>(Pcam, Vector3<Vector3<T>>(x0,x1,x2),Vector3<Vector3<T>>(x0c,x1c,x2c));
    }
    bool first_special_case1=true;
    template<class T> Data<T>
    special_case1(){
        // 2d points are orthogonal
return Data<T>();
        //turn Data<T>(Pcam, Vector3<Vector3<T>>(x0,x1,x2),Vector3<Vector3<T>>(x0c,x1c,x2c));
    }

};

}
