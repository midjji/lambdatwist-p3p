#include <utils/random.h>
namespace mlib{

/**
 * @brief randu integer random value drawn from uniform distribution
 * @param low  - included
 * @param high - included
 * @return random value according to compile time random value strategy
 */
double randu(double low, double high){
    std::uniform_real_distribution<double> rn(low,high);
    return rn(random::generator);
}

/**
 * @brief randn random value drawn from normal distribution
 * @param m
 * @param sigma
 * @return random value drawn from normal distribution
 */
double randn(double mean, double sigma){
    std::normal_distribution<double> rn(mean, sigma);
    return rn(random::generator);
}
/**
 * @brief randui integer random value drawn from uniform distribution
 * @param low  - included
 * @param high - included
 * @return random value according to compile time random value strategy
 */
int randui(int low, int high){
    std::uniform_int_distribution<int> rn(low,high);
    return rn(random::generator);
}




}
namespace cvl {




/// Uniform distribution on the unit sphere
template<class T,int R> cvl::Vector<T,R> getRandomUnitVector(){
    cvl::Vector<T,R> n=cvl::Vector<T,R>::Zero();
    // eliminating the inner sphere region has no effect on the manifold projection except to make it always valid
    while(n.squaredNorm() <1e-5)
        for(int i =0;i<R;++i)        n[i]=mlib::randn(0,1);

    n.normalize();
    return n;
}
// based on the above && the quaternion properties this gives a uniform distribution of rotations
Matrix3d random_rotation_matrix(){
    return getRotationMatrix(getRandomUnitVector<double,4>());
}
PoseD random_pose(){
    return PoseD(getRandomUnitVector<double,4>(), getRandomUnitVector<double,3>());
}

}
