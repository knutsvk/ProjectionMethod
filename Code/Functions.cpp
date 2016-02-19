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
    double alpha = 1+4*dt/Re*N*N;
    double beta = -dt/Re*N*N;

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

void updateLoadU(VectorXd u, VectorXd v, int N, double dt, 
        double a, double Re, VectorXd &f_U)
{
    /* Update the load vector used in the equation for the 
     * intermediate x-velocity. The stencil is such that each
     * point requires information from its own velocity from 
     * North, East, South and West, and from the y-velocity
     * from NorthEast, SouthEast, SouthWest and NorthWest. 
     * Boundary conditions from viscous term are included. */

    // CFL-related constant
    double alpha = dt*N; 

    // Viscosity-related constant
    double beta = dt/Re*N*N; 

    // For clarity, initiate doubles for each relevant point
    double u_0, u_N, u_E, u_S, u_W, v_NE, v_SE, v_SW, v_NW;
    
    for(int i=0; i<N-1; i++)
    {
        for(int j=0; j<N; j++)
        {
            // Point to be evaluated 
            u_0 = u[i*N+j];

            // x-velocities are zero if on west border
            if(i==0) u_W = 0.0;
            else u_W = u[(i-1)*N+j];

            // x-velocities are zero if on east border
            if(i==N-2) u_E = 0.0;
            else u_E = u[(i+1)*N+j];

            // y-velocities are zero if on south border
            // x-velocities from no-slip condition
            if(j==0)
            {
                u_S = -u_0;
                v_SE = 0.0;
                v_SW = 0.0;
            }
            else
            {
                u_S = u[i*N+(j-1)];
                v_SE = v[(j-1)*N+(i+1)];
                v_SW = v[(j-1)*N+i];
            }

            // y-velocities are zero if on north border
            // x-velocities from no-slip condition
            if(j==N-1)
            { 
                u_N = 2.0*a-u_0; 
                v_NE = 0.0;
                v_NW = 0.0;
            }
            else 
            {
                u_N = u[i*N+(j+1)];
                v_NE = v[j*N+(i+1)];
                v_NW = v[j*N+i];
            }
          
            // compute load vector
            f_U[i*N+j] = u_0-alpha*(
                    u_0*(u_E-u_W)+0.5*u_0*(v_NW+v_NE-v_SW-v_SE)
                    +0.125*(u_N-u_S)*(v_NW+v_NE+v_SE+v_SW));
            
            // add BCs from viscosiy term as appropriate
            // if on north or south boundary
            if(j==0)
                f_U[i*N+j] += beta*u_S;
            if(j==N-1)
                f_U[i*N+j] += beta*u_N;

        }
    }
}

void updateLoadV(VectorXd u, VectorXd v, int N, double dt, 
        double Re, VectorXd &f_V)
{
    /* Update the load vector used in the equation for the 
     * intermediate y-velocity. The stencil is such that each
     * point requires information from its own velocity from 
     * North, East, South and West, and from the x-velocity
     * from NorthEast, SouthEast, SouthWest and NorthWest. 
     * Boundary conditions from viscous term are included. */

    // CFL-related constant
    double alpha = dt*N; 

    // Viscosity-related constant
    double beta = dt/Re*N*N; 

    // For clarity, initiate doubles for each relevant point
    double v_0, v_N, v_E, v_S, v_W, u_NE, u_SE, u_SW, u_NW;
    
    for(int i=0; i<N-1; i++)
    {
        for(int j=0; j<N; j++)
        {
            // Point to be evaluated 
            v_0 = v[i*N+j];

            // y-velocities are zero if on south border
            if(i==0) v_S = 0.0;
            else v_S = v[(i-1)*N+j];

            // y-velocities are zero if on north border
            if(i==N-2) v_N = 0.0;
            else v_N = v[(i+1)*N+j];

            // x-velocities are zero if on west border
            // y-velocities from no-slip condition
            if(j==0)
            {
                v_W = -v_0;
                u_NW = 0.0;
                u_SW = 0.0;
            }
            else
            {
                v_W = v[i*N+(j-1)];
                u_NW = u[(j-1)*N+(i+1)];
                u_SW = u[(j-1)*N+i];
            }

            // x-velocities are zero if on east border
            // y-velocities from no-slip condition
            if(j==N-1)
            { 
                v_E = -v_0; 
                u_NE = 0.0;
                u_SE = 0.0;
            }
            else 
            {
                v_E = v[i*N+(j+1)];
                u_NE = u[j*N+(i+1)];
                u_SE = u[j*N+i];
            }
          
            // compute load vector
            f_V[i*N+j] = v_0-alpha*(
                    v_0*(v_N-v_S)+0.5*v_0*(u_NE+u_SE-u_NW-u_SW)
                    +0.125*(v_E-v_W)*(u_NW+u_NE+u_SE+u_SW));

            // add BCs from viscosiy term as appropriate
            // if on east or west boundary
            if(j==0)
                f_V[i*N+j] += beta*v_W;
            if(j==N-1)
                f_V[i*N+j] += beta*v_E;
        }
    }
}

void updateLoadp(VectorXd U, VectorXd V, int N, double dt, 
        VectorXd &f_p)
{
    /* Update the load vector used in the equation for the 
     * pressure. The stencil is such that each
     * point requires information from the intermediate 
     * x-velocities at
     * East and West and y-velocities at North and South. */

    // For clarity, initiate doubles for each relevant point
    double U_E, U_W, V_N, V_S;
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            // x-velocities are zero if on west border
            if(i==0) U_W = 0.0;
            else U_W = U[(i-1)*N+j];

            // x-velocities are zero if on east border
            if(i==N-1) U_E = 0.0;
            else U_E = U[i*N+j];

            // y-velocities are zero if on south border
            if(j==0) V_S = 0.0;
            else V_S = V[(j-1)*N+i];

            // y-velocities are zero if on north border
            if(j==N-1) V_N = 0.0;
            else V_N = V[j*N+i];
          
            // compute element of load vector
            f_p[i*N+j] = U_E-U_W+V_N-V_S;
        }
    }

    // Multiply by 1/CFL
    f_p *= 1.0/(N*dt);
}

void updateVelocities(VectorXd U, VectorXd V, VectorXd p, 
        int N, double dt, VectorXd &u, VectorXd &v)
{
    /* Update velocities after solving all three systems of 
     * equations. */

    // Approximation to the x-derivative of pressure
    VectorXd pDiffX = p.tail(N*(N-1))-p.head(N*(N-1));

    // Approximation to the y-deriavtive of pressure
    VectorXd pDiffY(N*(N-1));
/*    for(int i=0; i<N-1; i++)
    {
        pDiffY.segment(i*N, N) = 
            p.segment((i+1)*N,N)-p.segment(i*N,N);
    }*/

    for(int i=0; i<N-1; i++)
    {
        for(int j=0; j<N; j++)
        {
            pDiffY[i*N+j] = p[j*N+(i+1)]-p[j*N+i];
        }
    }

    // Compute
    u = U - dt*N*pDiffX;
    v = V - dt*N*pDiffY;
}

#endif
