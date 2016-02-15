#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <Eigen/Eigen>

using namespace std; 
using namespace Eigen;

MatrixXd buildVelocityMatrix(int N, double dt, double Re);
MatrixXd buildPressureMatrix(int N);

#endif
