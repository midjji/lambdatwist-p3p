#pragma once
#ifndef P3P_P3P_H
#define P3P_P3P_H
namespace kes{
class P3p {
 private:
  inline void polishQuarticRoots(const double *coeffs, double *roots);
  inline void vect_cross(const double *a, const double *b, double *result);
  inline double vect_dot(const double *a, const double *b);
  inline double vect_norm(const double *a);
  inline void vect_scale(const double s, const double *a, double *result);
  inline void vect_sub(const double *a, const double *b, double *result);
  inline void vect_divide(const double *a, const double d, double *result);
  inline void mat_mult(const double a[3][3], const double b[3][3], double result[3][3]);

 public:
  P3p() {}
  int computePoses(double **featureVectors, double **worldPoints, double solutions[3][16]);
  int solveQuartic(const double *factors, double *realRoots);
};
}//end namespace kes{

#endif //P3P_P3P_H
