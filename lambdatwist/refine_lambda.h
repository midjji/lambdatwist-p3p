#pragma once
#include <iostream>

namespace cvl{
template<class T>
void newton_refineL(Vector3<T>& L,
                    T a12, T a13, T a23,
                    T b12, T b13, T b23 ){
    for(int i=0;i<2;++i){
        T l1=L(0);
        T l2=L(1);
        T l3=L(2);
        T r1=l1*l1 + l2*l2 +b12*l1*l2 -a12;
        T r2=l1*l1 + l3*l3 +b13*l1*l3 -a13;
        T r3=l2*l2 + l3*l3 +b23*l2*l3 -a23;

        T dr1dl1=T(2.0)*l1 +b12*l2;
        T dr1dl2=T(2.0)*l2 +b12*l1;
        T dr1dl3=0;

        T dr2dl1=T(2.0)*l1 +b13*l3;
        T dr2dl2=0;
        T dr2dl3=T(2.0)*l3 +b13*l1;


        T dr3dl1=0;
        T dr3dl2=T(2.0)*l2 + b23*l3;
        T dr3dl3=T(2.0)*l3 + b23*l2;

        Vector3<T> g(r1*dr1dl1 + r2*dr2dl1 + r3*dr3dl1,
                     r1*dr1dl2 + r2*dr2dl2 + r3*dr3dl2,
                     r1*dr1dl3 + r2*dr2dl3 + r3*dr3dl3);
        Vector3<T> drdl(dr1dl1 + dr2dl1 +dr3dl1,
                        dr1dl2 + dr2dl2 +dr3dl2,
                        dr1dl3 + dr2dl3 +dr3dl3);
        Matrix<T,3,3> H=Matrix<T,3,3>(1,0,0,0,1,0,0,0,1)*1e-8 + drdl*drdl.transpose() + Matrix<T,3,3>(2.0*r1 + 2.0*r2, b12*r1,    b13*r2,
                                                                                                      b12*r1,          2*r1+2*r3, b23*r3,
                                                                                                      b13*r2,          b23*r3,    2*r2 +2*r3);
        L-=(H.inverse()*g);

    }

}
template<class T>
void sep_newton_refineL(Vector3<T>& L,
                        T a12, T a13, T a23,
                        T b12, T b13, T b23 ){

    //cout<<endl;
    //cout<<endl;

    //cout<<L<<endl;
    // l1
    for(int i=0;i<5;++i){
        {
            T l1=L(0);
            T l2=L(1);
            T l3=L(2);
            T r1=l1*l1 + l2*l2 +b12*l1*l2 -a12;
            T r2=l1*l1 + l3*l3 +b13*l1*l3 -a13;

            //cout<<r1<<" "<<r2<<endl;


            T dr1dl1=T(2.0)*l1 +b12*l2;
            T dr2dl1=T(2.0)*l1 +b13*l3;


            T dr1dl1dl1=2.0;
            T dr2dl1dl1=2.0;



            T dedl1=(dr1dl1*r1+dr2dl1*r2);
            T dedl1dl1=dr1dl1dl1*r1+dr2dl1dl1*r2 + dr1dl1*dr1dl1 + dr2dl1*dr2dl1;

            // should reduce the error!
            T delta=dedl1/dedl1dl1;
            //cout<<"delta: "<<delta<<endl;

            L(0)-= delta;


        }
        // l2
        {
            T l1=L(0);
            T l2=L(1);
            T l3=L(2);
            T r1=l1*l1 + l2*l2 +b12*l1*l2 -a12;
            T r3=l2*l2 + l3*l3 +b23*l2*l3 -a23;
            //cout<<r1<<" "<<r2<<endl;



            T dr1dl2=T(2.0)*l2 +b12*l1;
            T dr3dl2=T(2.0)*l2 + b23*l3;


            T dedl2=(dr1dl2*r1+dr3dl2*r3);
            T dedl2dl2=2*r1 +2*r3 + dr1dl2*dr1dl2 + dr3dl2*dr3dl2;

            // should reduce the error!
            T delta=dedl2/dedl2dl2;
            //cout<<"delta: "<<delta<<endl;

            L(1)-= delta;

        }
        // l3
        {
            T l1=L(0);
            T l2=L(1);
            T l3=L(2);
            T r2=l1*l1 + l3*l3 +b13*l1*l3 -a13;
            T r3=l2*l2 + l3*l3 +b23*l2*l3 -a23;

            T dr2dl3=T(2.0)*l3 +b13*l1;
            T dr3dl3=T(2.0)*l3 + b23*l2;


            T dedl3=(dr2dl3*r2+dr3dl3*r3);
            T dedl3dl3=2*r2 +2*r3 + dr2dl3*dr2dl3 + dr3dl3*dr3dl3;

            // should reduce the error!
            T delta=dedl3/dedl3dl3;
            //cout<<"delta: "<<delta<<endl;

            L(2)-= delta;

        }
    }
}


template<class T> Vector3<T> solve(const Matrix<T,3,3>& M,const  Vector3<T>& b){
    assert(M.is_symmetric());
    Vector3<T> x(0,0,0);
    Vector3<T> q;
    Vector3<T> r=b-M*x;
    T e=r.squaredLength();
    T e0,a,beta;
    Vector3<T> d=r;



    for(int i=0;i<3;++i){
        q=M*d;
        a=e/cvl::dot(q,d);
        x=x+a*d;
        r=r-a*q;
        e0=e;
        e=r.squaredLength();
        beta=e0/e;
        d=r+beta*d;
    }
    return x;
}




static int badsteps=0;
template<class T, int iterations>
/**
 * @brief refineL
 * @param L
 * @param a12
 * @param a13
 * @param a23
 * @param b12
 * @param b13
 * @param b23
 *
 * Gauss-Newton Solver
 * For unknown reasons it always works for the correct solution, but not always for the other solutions!
 *
 */
void gauss_newton_refineL(Vector3<T>& L,
                          T a12, T a13, T a23,
                          T b12, T b13, T b23 ){

    // const expr makes it easier for the compiler to unroll
    for(int i=0;i<iterations;++i){
        T l1=L(0);
        T l2=L(1);
        T l3=L(2);
        T r1=l1*l1 + l2*l2 +b12*l1*l2 -a12;
        T r2=l1*l1 + l3*l3 +b13*l1*l3 -a13;
        T r3=l2*l2 + l3*l3 +b23*l2*l3 -a23;

        if(std::abs(r1) +std::abs(r2) +std::abs(r3)<1e-10) break;




        T dr1dl1=(2.0)*l1 +b12*l2;
        T dr1dl2=(2.0)*l2 +b12*l1;
        //T dr1dl3=0;

        T dr2dl1=(2.0)*l1 +b13*l3;
        //T dr2dl2=0;
        T dr2dl3=(2.0)*l3 +b13*l1;


        //T dr3dl1=0;
        T dr3dl2=(2.0)*l2 + b23*l3;
        T dr3dl3=(2.0)*l3 + b23*l2;



        Vector3<T> r(r1, r2, r3);
        /*


*/
        // or skip the inverse and make it explicit...
        {
#if 1
            T v0=dr1dl1;
            T v1=dr1dl2;
            T v3=dr2dl1;
            T v5=dr2dl3;
            T v7=dr3dl2;
            T v8=dr3dl3;
            T det=(1.0)/(- v0*v5*v7 - v1*v3*v8);

            Matrix<T,3,3> Ji( -v5*v7, -v1*v8,  v1*v5,
                              -v3*v8,  v0*v8, -v0*v5,
                              v3*v7, -v0*v7, -v1*v3);
            Vector3<T> L1=Vector3<T>(L) - det*(Ji*r);
#else

            Matrix<T,3,3> J(dr1dl1, dr1dl2, 0,
                            dr2dl1,0, dr2dl3,
                            0, dr3dl2, dr3dl3);

            //%l=l - g*H\G;%inv(H)*G
            //L=L - g*J\r; //% works because the size is ok!
#if 1
            Vector3<T> L1=Vector3<T>(L) - (J.inverse()*r);

#else
            Vector3<T> L1=Vector3<T>(L) - solve(J.transpose()*J,J.transpose()*r);
#endif

#endif


            {



                T l1=L(0);
                T l2=L(1);
                T l3=L(2);
                T r11=l1*l1 + l2*l2 +b12*l1*l2 -a12;
                T r12=l1*l1 + l3*l3 +b13*l1*l3 -a13;
                T r13=l2*l2 + l3*l3 +b23*l2*l3 -a23;
                if(std::abs(r11) +std::abs(r12) + std::abs(r13)>std::abs(r1) +std::abs(r2) +std::abs(r3)){
                    // cout<<"bad step: "<< det*(Ji*r)<<++badsteps<<" "<< i<<endl;
                    break;
                }
                else
                    L=L1;
            }
        }




    }
    // cout<<i<<endl;


}

template<class T> inline Vector3<T> get_error(Vector3<T>& L,
                                       T a12, T a13, T a23,
                                       T b12, T b13, T b23){
    T l1=L(0);
    T l2=L(1);
    T l3=L(2);
    T r1=l1*l1 + l2*l2 +b12*l1*l2 -a12;
    T r2=l1*l1 + l3*l3 +b13*l1*l3 -a13;
    T r3=l2*l2 + l3*l3 +b23*l2*l3 -a23;

    return Vector3<T> (r1, r2, r3);

}

template<class T> Matrix<T,3,3> get_J(Vector3<T>& L,
                                      T a12, T a13, T a23,
                                      T b12, T b13, T b23){
    T l1=L(0);
    T l2=L(1);
    T l3=L(2);
    T dr1dl1=(2.0)*l1 +b12*l2;
    T dr1dl2=(2.0)*l2 +b12*l1;
    //T dr1dl3=0;

    T dr2dl1=(2.0)*l1 +b13*l3;
    //T dr2dl2=0;
    T dr2dl3=(2.0)*l3 +b13*l1;


    //T dr3dl1=0;
    T dr3dl2=(2.0)*l2 + b23*l3;
    T dr3dl3=(2.0)*l3 + b23*l2;
    return Matrix<T,3,3>(dr1dl1, dr1dl2, 0,
                         dr2dl1,0, dr2dl3,
                         0, dr3dl2, dr3dl3);



}
template<class T> inline Vector3<T> get_delta(Vector3<T>& L,
                            T a12, T a13, T a23,
                            T b12, T b13, T b23,
                                       Vector3<T>& r){
    T l1=L(0);
    T l2=L(1);
    T l3=L(2);
    T dr1dl1=(2.0)*l1 +b12*l2;
    T dr1dl2=(2.0)*l2 +b12*l1;
    //T dr1dl3=0;

    T dr2dl1=(2.0)*l1 +b13*l3;
    //T dr2dl2=0;
    T dr2dl3=(2.0)*l3 +b13*l1;


    //T dr3dl1=0;
    T dr3dl2=(2.0)*l2 + b23*l3;
    T dr3dl3=(2.0)*l3 + b23*l2;
    T v0=dr1dl1;
    T v1=dr1dl2;
    T v3=dr2dl1;
    T v5=dr2dl3;
    T v7=dr3dl2;
    T v8=dr3dl3;
    T det=(1.0)/(- v0*v5*v7 - v1*v3*v8);

    Matrix<T,3,3> Ji( -v5*v7, -v1*v8,  v1*v5,
                      -v3*v8,  v0*v8, -v0*v5,
                      v3*v7, -v0*v7, -v1*v3);
    return - det*(Ji*r);
}


template<class T, int iterations>
/**
 * @brief refineL
 * @param L
 * @param a12
 * @param a13
 * @param a23
 * @param b12
 * @param b13
 * @param b23
 *
 * Gauss-Newton Solver
 * For unknown reasons it always works for the correct solution, but not always for the other solutions!
 *
 */
void gauss_newton_refineL2(Vector3<T>& L,
                           T a12, T a13, T a23,
                           T b12, T b13, T b23 ){

    //slightly worse than the basic one

    // const expr makes it easier for the compiler to unroll
    Vector3<T> r1,L1;
    Vector3<T> r(get_error(L,a12,a13,a23,b12,b13,b23));
    T e=r.absSum();
    T e1=r.absSum()+1;
    T a=1;int inner=0;
    for(int i=0;i<iterations;++i){


        if(e<1e-10) break;

        Vector3<T> delta(get_delta(L,a12,a13,a23,b12,b13,b23,r));
        while(e1>e){
            L1=(L + a*delta);
            r1=get_error(L1,a12,a13,a23,b12,b13,b23);
            T e1=r1.absSum();
            if(e1>e)
                a=0.5*a;
            inner++;
            if(inner>5)
                break;
        }
        a=1;
        inner=0;
        r=r1;
        L=L1;
        e=e1;
    }
    // cout<<i<<endl;


}

}
