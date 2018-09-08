#pragma once
#include <iostream>
#include <iomanip>


#include <Eigen/Dense>
#include <utils/cvl/pose.h>
#include <utils/matlab_helpers.h>

namespace cvl {


template<class T>
class  Data{
public:
    cvl::Vector3<Vector3<T>>  x0,xr;
    Pose<T> P;

    Data(){}

    Data(Pose<T> p,
         Vector3<Vector3<T> > x0,
         Vector3<Vector3<T> > xr):x0(x0),xr(xr),P(p){}


    Vector3<T> getLambda()
    {
        return {xr[0].length(),xr[1].length(),xr[2].length()};
    }

    cvl::Vector3<Vector3<T> > getYh() const
    {
        Vector3<Vector3<T>> y;
        for(int i=0;i<3;++i)
            y[i]=xr[i].normalized();

        return y;
    }

    Eigen::Matrix3d getFeatureVectors() const {
        Eigen::Matrix3d feature_vectors;
        // should possibly be normalized...?
        Vector3<Vector3<T>> ys;
        for(int i=0;i<3;++i)
            ys[i]=xr[i].normalized();
        feature_vectors<<
                          ys[0][0],ys[1][0],ys[2][0],
                ys[0][1],ys[1][1],ys[2][1],
                ys[0][2],ys[1][2],ys[2][2];

        return feature_vectors;
    }



    Eigen::Matrix3d getWorldPoints() const {
        Eigen::Matrix3d world_points;


        world_points<<
                       x0[0][0],x0[1][0],x0[2][0],
                x0[0][1],x0[1][1],x0[2][1],
                x0[0][2],x0[1][2],x0[2][2];
        return world_points;
    }
    std::string toMatlab() const {
        std::stringstream ss;

        ss<<std::setprecision(20)<<"data.pose="<<mlib::getMatlabMatrix(P,20)<<"\n";
        Vector3<Vector3<T>> ys=getYh();
        for(int j=0;j<3;++j)
            ss<<"data.y"<<j<<"=["<<ys[j][0]<<" "<<ys[j][1]<<" "<<ys[j][2]<<"]';\n";
        for(int j=0;j<3;++j)
            ss<<"data.x"<<j<<"=["<<x0[j][0]<<" "<<x0[j][1]<<" "<<x0[j][2]<<"]';\n";
        for(int j=0;j<3;++j)
            ss<<"data.xr"<<j<<"=["<<xr[j][0]<<" "<<xr[j][1]<<" "<<xr[j][2]<<"]';\n";
        return ss.str();
    }

    T min_error(cvl::Vector<cvl::Matrix<T,3,3>,4>& Rs,
                cvl::Vector<cvl::Vector3<T>,4>& Ts, int valid){

        T min_error=1; // 1 is effectively a failure anyways, except for long Ts
        for(int i=0;i<valid;++i){
            // the square is effectively the max norm otherwize for small values...
            T error=(Rs[i] - P.getR()).abs().sum() + (Ts[i] -P.getT()).abs().sum();
            if(!std::isnan(error))
                min_error=std::min(error,min_error);
        }
        return min_error;
    }

    bool good_solutions(cvl::Matrix<T,3,3> R,
                        cvl::Vector3<T> t, bool full=false)  const{
        // note there may be more valid solutions beyond the one which generated the sample!

        // verify not nan, kneip fails on this one...
        //if((R.isnan()  || t.isnan())) return false;
        // verify points infront , kneip and ke fail this one
        for(int j=0;j<3;++j)             if((R*x0[j] +t)[2]<0) return false;

        if(full){
            // verify rotation matrix
            if(std::abs(R.determinant()-1)>1e-6) return false;
            // verify rotation matrix
            if(((R.transpose()*R) -cvl::Matrix<T,3,3>(1,0,0,0,1,0,0,0,1)  ).abs().sum()>1e-6) return false;

            // verify the repojection error
            T err=0;
            for(int j=0;j<3;++j)  err+=((R*x0[j] +t).dehom() - (xr[j].dehom())).abs().sum();
            if(err>1e-5) return false;
        }
        return true;

    }



    int good_solutions(cvl::Vector<cvl::Matrix<T,3,3>,4>& Rs,
                       cvl::Vector<cvl::Vector3<T>,4>& Ts, int valid, int& duplicates){
        // all valid solutions must be in the start of the Vector, no invalid solutions after the first bad one!

        int goods=0;

        std::vector<int> valids;
        for(int i=0;i<valid;++i){
            if(!good_solutions(Rs[i],Ts[i],true)) continue;

            goods++;
            bool dup=false;

            for(int j:valids){

                if((Rs[i]-Rs[j]).abs().sum() + (Ts[i]-Ts[j]).abs().sum()<1e-5){
                    dup=true;
                }
            }

            if(!dup)
                valids.push_back(i);
            else
                duplicates++;
        }


        return goods;
    }

};



}
