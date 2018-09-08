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

class Generator{
public:

    Generator(){

    }
    template<class T>
    Data<T> next(){

        cvl::Matrix<T,3,3> R=mlib::getRandomRotation<double>();
        //cout<<"Warning should it be with a mean of -1?"<<endl;
        //cvl::Vector3<T> t(mlib::randn<long double>(-1,1),mlib::randn<long double>(-1,1),mlib::randn<long double>(-1,1));
        // slower to generate...
        cvl::Vector3<T> t(mlib::randn<long double>(0,1),mlib::randn<long double>(0,1),mlib::randn<long double>(0,1));
        //std::cout<<"Data:\n"<<R<<" "<<t<<std::endl;
          //t*=mlib::randu<long double>(0.1,1000);// doesnt matter
        //t*=0.01;
        Pose<T> P(R,t);
        Vector3<Vector3<T>>  x0,xr;
        Vector3<Vector2<T>>  ys;



        // good generator, pick a uniform position between -1,1 and multiply by a distance...

        for(int i=0;i<3;++i){
            //Vector2d y((double)((float)randn(0,1)),(double)((float)randn(0,1)));
            Vector2<T> y(mlib::randu<long double>(-1,1),mlib::randu<long double>(-1,1));
          // y*=0.1;
            T z=mlib::randu<long double>(0.1,100);// 0.1 to 10 is good
            //T z=mlib::randu<long double>(0.9,1.1);// 0.1 to 1000 is good too

            Vector3<T> xri=z*y.homogeneous();
            Vector3<T> x=(P.inverse()*xri);
            x0[i]=x;
            xr[i]=xri;

           // y+=Vector2<T>(mlib::randn<long double>(0.0,0.01),mlib::randn<long double>(0,0.01));
            ys[i]=y;
        }


        // the reconstructed points must be inside the -1,1 of the second image...
        bool visible_by_both=true;
        for(auto x:x0){
            auto y=x.dehom();
            if(y.absMax()>1) visible_by_both=false;
        }
        visible_by_both=true; // there isnt really a second image... still sometimes some x might wind up beeing gigantic and that is not intentional...
        // they must not be exactly colinear
        bool colinear=false;
        {
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    if(i!=j){
                        if(x0(i)==x0(j) || ys[i]==ys[j])
                            colinear=true;
                    }
        }
        if(!visible_by_both || colinear) return next<T>();





        // minimum dist between obs?
        //bool dist=(ys[0]-ys[1]).length()<1e-2 || (ys[0]-ys[2]).length()<1e-2 || (ys[1]-ys[2]).length()<1e-2;
        //bool dist2=(x0[0]-x0[1]).length()<1e-2 || (x0[0]-x0[2]).length()<1e-2 || (x0[1]-x0[2]).length()<1e-2;


        // minimum angle between the obs
        //if(abs(ys[0].dot(ys[1]))<1e-4 || abs(ys[0].dot(ys[2]))<1e-4||abs(ys[2].dot(ys[1]))<1e-4)
        //  dist2=true;
        // minimum angle between the coordinates

        Vector3<T> a=xr[1]-xr[0];
        a.normalize();
        Vector3<T> b=xr[2]-xr[1];
        b.normalize();
        bool dist4=a.cross(b).length()<1e-3;
        a=ys[0].homogeneous()-ys[1].homogeneous();
        b=ys[0].homogeneous()-ys[2].homogeneous();
        a.normalize();
        b.normalize();
        dist4=dist4 || a.cross(b).length()<1e-3;

        //std::cout<<dist<<" "<<dist2<< " "<<dist4<<std::endl;

        //if(dist4|| dist ||dist2)            return next<T>();


        // reorder them in a particular order? say by angle diff so the second one has max angle to the others...

        //if(angle(ys[0],ys[1])>)
        return Data<T>(P,x0,xr);
    }



};

}
