#include "kneip/kneip.h"


int P3P::computePoses(const Eigen::Matrix3d & feature_vectors, const Eigen::Matrix3d & world_points,
                      Eigen::Matrix<Eigen::Matrix<double, 3, 4>, 4, 1> & solutions)
{
  // Extraction of world points
  Eigen::Vector3d P1 = world_points.col(0);
  Eigen::Vector3d P2 = world_points.col(1);
  Eigen::Vector3d P3 = world_points.col(2);

  // Verification that world points are not colinear
  Eigen::Vector3d temp1 = P2 - P1;
  Eigen::Vector3d temp2 = P3 - P1;

  if (temp1.cross(temp2).norm() == 0)  {    return -1;  }

  // Extraction of feature vectors
  Eigen::Vector3d f1 = feature_vectors.col(0);
  Eigen::Vector3d f2 = feature_vectors.col(1);
  Eigen::Vector3d f3 = feature_vectors.col(2);

  // Creation of intermediate camera frame
  Eigen::Vector3d e1 = f1;
  Eigen::Vector3d e3 = f1.cross(f2);
  e3 = e3 / e3.norm();
  Eigen::Vector3d e2 = e3.cross(e1);

  Eigen::Matrix3d T;
  T.row(0) = e1.transpose();
  T.row(1) = e2.transpose();
  T.row(2) = e3.transpose();

  f3 = T * f3;

  // Reinforce that f3(2,0) > 0 for having theta in [0;pi]
  if (f3(2, 0) > 0)
  {
    f1 = feature_vectors.col(1);
    f2 = feature_vectors.col(0);
    f3 = feature_vectors.col(2);

    e1 = f1;
    e3 = f1.cross(f2);
    e3 = e3 / e3.norm();
    e2 = e3.cross(e1);

    T.row(0) = e1.transpose();
    T.row(1) = e2.transpose();
    T.row(2) = e3.transpose();

    f3 = T * f3;

    P1 = world_points.col(1);
    P2 = world_points.col(0);
    P3 = world_points.col(2);
  }

  // Creation of intermediate world frame
  Eigen::Vector3d n1 = P2 - P1;
  n1 = n1 / n1.norm();
  Eigen::Vector3d n3 = n1.cross(P3 - P1);
  n3 = n3 / n3.norm();
  Eigen::Vector3d n2 = n3.cross(n1);

  Eigen::Matrix3d N;
  N.row(0) = n1.transpose();
  N.row(1) = n2.transpose();
  N.row(2) = n3.transpose();

  // Extraction of known parameters
  P3 = N * (P3 - P1);

  double d_12 = (P2 - P1).norm();
  double f_1 = f3(0, 0) / f3(2, 0);
  double f_2 = f3(1, 0) / f3(2, 0);
  double p_1 = P3(0, 0);
  double p_2 = P3(1, 0);

  double cos_beta = f1.dot(f2);
  double b = 1 / (1 - pow(cos_beta, 2)) - 1;

  if (cos_beta < 0)
  {
    b = -sqrt(b);
  }
  else
  {
    b = sqrt(b);
  }

  // Definition of temporary variables for avoiding multiple computation
  double f_1_pw2 = pow(f_1, 2);
  double f_2_pw2 = pow(f_2, 2);
  double p_1_pw2 = pow(p_1, 2);
  double p_1_pw3 = p_1_pw2 * p_1;
  double p_1_pw4 = p_1_pw3 * p_1;
  double p_2_pw2 = pow(p_2, 2);
  double p_2_pw3 = p_2_pw2 * p_2;
  double p_2_pw4 = p_2_pw3 * p_2;
  double d_12_pw2 = pow(d_12, 2);
  double b_pw2 = pow(b, 2);

  // Computation of factors of 4th degree polynomial
  Eigen::Matrix<double, 5, 1> factors;

  factors(0) = -f_2_pw2 * p_2_pw4 - p_2_pw4 * f_1_pw2 - p_2_pw4;

  factors(1) = 2 * p_2_pw3 * d_12 * b + 2 * f_2_pw2 * p_2_pw3 * d_12 * b - 2 * f_2 * p_2_pw3 * f_1 * d_12;

  factors(2) = -f_2_pw2 * p_2_pw2 * p_1_pw2 - f_2_pw2 * p_2_pw2 * d_12_pw2 * b_pw2 - f_2_pw2 * p_2_pw2 * d_12_pw2
      + f_2_pw2 * p_2_pw4 + p_2_pw4 * f_1_pw2 + 2 * p_1 * p_2_pw2 * d_12 + 2 * f_1 * f_2 * p_1 * p_2_pw2 * d_12 * b
      - p_2_pw2 * p_1_pw2 * f_1_pw2 + 2 * p_1 * p_2_pw2 * f_2_pw2 * d_12 - p_2_pw2 * d_12_pw2 * b_pw2
      - 2 * p_1_pw2 * p_2_pw2;

  factors(3) = 2 * p_1_pw2 * p_2 * d_12 * b + 2 * f_2 * p_2_pw3 * f_1 * d_12 - 2 * f_2_pw2 * p_2_pw3 * d_12 * b
      - 2 * p_1 * p_2 * d_12_pw2 * b;

  factors(4) = -2 * f_2 * p_2_pw2 * f_1 * p_1 * d_12 * b + f_2_pw2 * p_2_pw2 * d_12_pw2 + 2 * p_1_pw3 * d_12
      - p_1_pw2 * d_12_pw2 + f_2_pw2 * p_2_pw2 * p_1_pw2 - p_1_pw4 - 2 * f_2_pw2 * p_2_pw2 * p_1 * d_12
      + p_2_pw2 * f_1_pw2 * p_1_pw2 + f_2_pw2 * p_2_pw2 * d_12_pw2 * b_pw2;

  // Computation of roots
  Eigen::Matrix<double, 4, 1> realRoots;

  P3P::solveQuartic(factors, realRoots);

  // Backsubstitution of each solution
  for (int i = 0; i < 4; ++i)
  {
    double cot_alpha = (-f_1 * p_1 / f_2 - realRoots(i) * p_2 + d_12 * b)
        / (-f_1 * realRoots(i) * p_2 / f_2 + p_1 - d_12);

    double cos_theta = realRoots(i);
    double sin_theta = sqrt(1 - pow((double)realRoots(i), 2));
    double sin_alpha = sqrt(1 / (pow(cot_alpha, 2) + 1));
    double cos_alpha = sqrt(1 - pow(sin_alpha, 2));

    if (cot_alpha < 0)
    {
      cos_alpha = -cos_alpha;
    }

    Eigen::Vector3d C;
    C(0) = d_12 * cos_alpha * (sin_alpha * b + cos_alpha);
    C(1) = cos_theta * d_12 * sin_alpha * (sin_alpha * b + cos_alpha);
    C(2) = sin_theta * d_12 * sin_alpha * (sin_alpha * b + cos_alpha);

    C = P1 + N.transpose() * C;

    Eigen::Matrix3d R;
    R(0, 0) = -cos_alpha;
    R(0, 1) = -sin_alpha * cos_theta;
    R(0, 2) = -sin_alpha * sin_theta;
    R(1, 0) = sin_alpha;
    R(1, 1) = -cos_alpha * cos_theta;
    R(1, 2) = -cos_alpha * sin_theta;
    R(2, 0) = 0;
    R(2, 1) = -sin_theta;
    R(2, 2) = cos_theta;

    R = N.transpose() * R.transpose() * T;

    Eigen::Matrix<double, 3, 4> solution;
    solution.block<3, 3>(0, 0) = R;
    solution.col(3) = C;

    solutions(i) = solution;
  }

  return 0;
}

int P3P::solveQuartic(const Eigen::Matrix<double, 5, 1> & factors, Eigen::Matrix<double, 4, 1> & real_roots)
{
  double A = factors(0);
  double B = factors(1);
  double C = factors(2);
  double D = factors(3);
  double E = factors(4);

  double A_pw2 = A * A;
  double B_pw2 = B * B;
  double A_pw3 = A_pw2 * A;
  double B_pw3 = B_pw2 * B;
  double A_pw4 = A_pw3 * A;
  double B_pw4 = B_pw3 * B;

  double alpha = -3 * B_pw2 / (8 * A_pw2) + C / A;
  double beta = B_pw3 / (8 * A_pw3) - B * C / (2 * A_pw2) + D / A;
  double gamma = -3 * B_pw4 / (256 * A_pw4) + B_pw2 * C / (16 * A_pw3) - B * D / (4 * A_pw2) + E / A;

  double alpha_pw2 = alpha * alpha;
  double alpha_pw3 = alpha_pw2 * alpha;

  std::complex<double> P(-alpha_pw2 / 12 - gamma, 0);
  std::complex<double> Q(-alpha_pw3 / 108 + alpha * gamma / 3 - pow(beta, 2) / 8, 0);
  std::complex<double> R = -Q / 2.0 + sqrt(pow(Q, 2.0) / 4.0 + pow(P, 3.0) / 27.0);

  std::complex<double> U = pow(R, (1.0 / 3.0));
  std::complex<double> y;

  if (U.real() == 0)
    y = -5.0 * alpha / 6.0 - pow(Q, (1.0 / 3.0));
  else
    y = -5.0 * alpha / 6.0 - P / (3.0 * U) + U;

  std::complex<double> w = sqrt(alpha + 2.0 * y);

  std::complex<double> temp;

  temp = -B / (4.0 * A) + 0.5 * (w + sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
  real_roots(0) = temp.real();
  temp = -B / (4.0 * A) + 0.5 * (w - sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
  real_roots(1) = temp.real();
  temp = -B / (4.0 * A) + 0.5 * (-w + sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
  real_roots(2) = temp.real();
  temp = -B / (4.0 * A) + 0.5 * (-w - sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
  real_roots(3) = temp.real();

  return 0;
}


