#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std; 
using namespace Eigen;

MatrixXd buildVelocityMatrix(int N, double dt, double Re);
MatrixXd buildPressureMatrix(int N);
void updateLoadU(VectorXd u, VectorXd v, int N, double dt, 
        double a, double Re, VectorXd &f_U);
void updateLoadV(VectorXd u, VectorXd v, int N, double dt, 
        double a, double Re, VectorXd &f_V);
void updateLoadp(VectorXd U, VectorXd V, int N, double dt, 
        VectorXd &f_p);
void updateVelocities(VectorXd U, VectorXd V, VectorXd p, 
        int N, double dt, VectorXd &u, VectorXd &v);

#endif
