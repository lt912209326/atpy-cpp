
#ifndef _TWISS_H_
#define _TWISS_H_

#include <Eigenvalues>
#include <Dense>
#include <exception>
#include "cppconstants.h"
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::EigenSolver;
using Eigen::Map;
using Eigen::RowMajor;


typedef Eigen::Matrix<double,2,2,RowMajor > Matrix2d;
typedef Eigen::Matrix<double,4,4,RowMajor > Matrix4d;
typedef Eigen::Matrix<double,6,6,RowMajor > Matrix6d;

Matrix2d _symplectic_conjugate_2by2(Matrix2d& m);

int calculate_coupled_period_twiss(double* mat, double* tws);

int propagate_twiss(double* tout, double* tin, double* mat);

double  _propagate_decoupled_twiss_2x2(double* tout,double* tin, Matrix2d& mat);

inline double phase_cumsum(double nu1,double nu0){
    double frac_nu0=fmod(nu0,1.0);
    return frac_nu0<nu1? nu0-frac_nu0+nu1 : nu0-frac_nu0+nu1+1 ;
}

#endif
