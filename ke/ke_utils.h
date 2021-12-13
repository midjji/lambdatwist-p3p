#pragma once

#include <ke/ke.h>
#include <utils/cvl/matrix.h>
#include <vector>
#include <data.h>
#include <solver.h>
namespace kes{


template<class T>
int ke_p3p_fair(const cvl::Data<T>& data,
                cvl::Vector<cvl::Matrix<T,3,3>,4>& Rs,
                cvl::Vector<cvl::Vector3<T>,4>& Ts){



    // all of this datashuffling is optimized away at O3 and Ofast, check the asm if you think otherwize...
    cvl::Vector3d yh0=data.xr[0];
    cvl::Vector3d yh1=data.xr[1];
    cvl::Vector3d yh2=data.xr[2];
    yh0.normalize();
    yh1.normalize();
    yh2.normalize();

    //double y0[3]={yh0[0],yh0[1],yh0[2]};
    //double y1[3]={yh1[0],yh1[1],yh1[2]};
    //double y2[3]={yh2[0],yh2[1],yh2[2]};
    double y0[3]={yh0[0],yh1[0],yh2[0]};
    double y1[3]={yh0[1],yh1[1],yh2[1]};
    double y2[3]={yh0[2],yh1[2],yh2[2]};
    double* yy[3];
    yy[0]=&y0[0];
    yy[1]=&y1[0];
    yy[2]=&y2[0];

    cvl::Vector3d xw0=data.x0[0];
    cvl::Vector3d xw1=data.x0[1];
    cvl::Vector3d xw2=data.x0[2];

    double x0[3]={xw0[0],xw1[0],xw2[0]};
    double x1[3]={xw0[1],xw1[1],xw2[1]};
    double x2[3]={xw0[2],xw1[2],xw2[2]};
    double* xx[3];

    xx[0]=&x0[0];
    xx[1]=&x1[0];
    xx[2]=&x2[0];



    kes::P3p p3p;
    double solutions[3][16];
    p3p.computePoses(yy,xx,solutions);
    int sols=0;
    cvl::Matrix<T,3,3> R;
    cvl::Vector<T,3> t;
    // they need to be converted to a standard format... again, all is basically optimized away except the tests
    for(int i=0;i<4;++i){

        R(0,0)=solutions[0][i * 4 + 1];// = R[0][0];
        R(1,0)=solutions[1][i * 4 + 1];// = R[1][0];
        R(2,0)=solutions[2][i * 4 + 1];// = R[2][0];
        R(0,1)=solutions[0][i * 4 + 2];// = R[0][1];
        R(1,1)=solutions[1][i * 4 + 2];// = R[1][1];
        R(2,1)=solutions[2][i * 4 + 2];// = R[2][1];
        R(0,2)=solutions[0][i * 4 + 3];// = R[0][2];
        R(1,2)=solutions[1][i * 4 + 3];// = R[1][2];
        R(2,2)=solutions[2][i * 4 + 3];// = R[2][2];

        t[0] = solutions[0][i * 4];// = temp[0];
        t[1] = solutions[1][i * 4];// = temp[1];
        t[2] = solutions[2][i * 4];// = temp[2];


        // change to camera centered for testing
        Rs[sols]=R.transpose();
        Ts[sols]=-Rs[sols]*t;
        // this post-processing must always be performed since no one would want invalid or incorrect solutions.
        // the thresholds are generally rather generous and could easily be lower
        // total time is roughly 40-60ns per call to Ke's... or 140 on ubuntu 18?

        if(data.good_solutions(Rs[sols],Ts[sols])) sols++;

    }
    // exit(1);
    return sols;
}








}
