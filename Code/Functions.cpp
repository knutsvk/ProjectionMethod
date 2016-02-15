#ifndef __FUNCTIONS_CPP
#define __FUNCTIONS_CPP

#include "Functions.h"

MatrixXd buildVelocityMatrix(int N, double dt, double Re)
{
    /* Builds the stiffness matrix used to solve the system of
     * equations Ax=b when x is an intermediate velocity. 
     * A is a N*(N-1)-by-N*(N-1) block matrix, consisting of 
     * diagonal blocks B which themselves are tridiagonal and 
     * subdiagonal blocks C which are diagonal. */

    // Constants for the diagonals
    double alpha = 1+4*dt/Re;
    double beta = -dt/Re;

    // Build the tridiagonal matrix which constitutes a single
    // block on the diagonal of A 
    MatrixXd B = MatrixXd::Zero(N,N);
    B.diagonal() = alpha*VectorXd::Ones(N);
    B.diagonal(1) = beta*VectorXd::Ones(N-1);
    B.diagonal(-1) = beta*VectorXd::Ones(N-1);

    // Build the diagonal matrix which is on the subdiagonals
    MatrixXd C = MatrixXd::Zero(N,N);
    C.diagonal() = beta*VectorXd::Ones(N);

    // Populate A with the blocks
    MatrixXd A = MatrixXd::Zero(N*(N-1),N*(N-1));
    for(int i=0; i<N-1; i++)
    {
        A.block(i*N,i*N,N,N) = B;
        if(i>0)
            A.block((i-1)*N,i*N,N,N) = C;
        if(i<N-2)
            A.block((i+1)*N,i*N,N,N) = C;
    }

    return A;
}

MatrixXd buildPressureMatrix(int N)
{
    /* Builds the stiffness matrix used to solve the system of
     * equations Ax=b when x is the pressure. 
     * A is a N*N-by-N*N block matrix, consisting of 
     * diagonal blocks B which themselves are tridiagonal and 
     * subdiagonal blocks C which are diagonal. */

    // Build the tridiagonal matrix which constitutes a single
    // block on the diagonal of A 
    MatrixXd B = MatrixXd::Zero(N,N);
    B.diagonal() = -4*VectorXd::Ones(N);
    B(0,0) = -3;
    B(N-1,N-1) = -3;
    B.diagonal(1) = VectorXd::Ones(N-1);
    B.diagonal(-1) = VectorXd::Ones(N-1);

    // Build the diagonal matrix which is on the subdiagonals
    MatrixXd C = MatrixXd::Identity(N,N);

    // Populate A with the blocks
    MatrixXd A = MatrixXd::Zero(N*N,N*N);
    for(int i=0; i<N; i++)
    {
        A.block(i*N,i*N,N,N) = B;
        if(i==0)
            A.block(i*N,i*N,N,N) += C;
        else 
            A.block((i-1)*N,i*N,N,N) = C;
        if(i==N-1)
            A.block(i*N,i*N,N,N) += C;
        else
            A.block((i+1)*N,i*N,N,N) = C;
    }

    return A;
}

#endif
