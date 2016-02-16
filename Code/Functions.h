#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <Eigen/Eigen>

using namespace std; 
using namespace Eigen;

MatrixXd buildVelocityMatrix(int N, double dt, double Re);
MatrixXd buildPressureMatrix(int N);
void updateLoadU(VectorXd u, VectorXd v, int N, double dt, 
        double a, double Re, VectorXd &f_U);
void updateLoadV(VectorXd u, VectorXd v, int N, double dt, 
        double a, double Re, VectorXd &f_V);

#endif
