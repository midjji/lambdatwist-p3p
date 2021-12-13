#pragma once
#include <utils/cvl/matrix.h>
#include <vector>
#include <data.h>
namespace kneip{




template<class T>
int kneip_p3p_fair(const cvl::Data<T>& data,
                   cvl::Vector<cvl::Matrix<T,3,3>,4>& Rs,
                   cvl::Vector<cvl::Vector3<T>,4>& Ts){

    Eigen::Matrix<Eigen::Matrix<double, 3, 4>, 4, 1>  solutions;

    P3P::computePoses(data.getFeatureVectors(),data.getWorldPoints(),solutions);



    int sols=0;
    for(int i=0;i<4;++i){
        Eigen::Matrix<double,3,4> a=solutions(i);

        cvl::Matrix<T,3,3> R(a(0,0),a(0,1),a(0,2),
                             a(1,0),a(1,1),a(1,2),
                             a(2,0),a(2,1),a(2,2));
        cvl::Vector3<T> t(a(0,3),a(1,3),a(2,3));

        Rs[sols]=R.transpose();
        Ts[sols]=-Rs[sols]*t;

        if(data.good_solutions(Rs[sols],Ts[sols])) sols++;
    }
    return sols;
}

}// end namespace kneip
