#ifndef __FUNCTIONS_CPP
#define __FUNCTIONS_CPP

#include "Functions.h"

VectorXd doubleShearLayer(int N)
{
    double y; 
    double dx = 1.0 / N; 
    VectorXd u(N * N);
    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            y = ( 0.5 + j ) * dx; 
            if( y <= 0.5 ) u[i * N + j] = tanh( 30 * ( y - 0.25 ) );
            else u[i * N + j] = tanh( 30 * ( 0.75 - y ) );
        }
    }
    return u;
}

VectorXd smallPerturbation(int N)
{
    double x; 
    double dx = 1.0 / N; 
    VectorXd v(N * N);
    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            x = ( 0.5 + i ) * dx; 
            v[j * N + i] = 0.05 * sin( 2 * M_PI * x );
        }
    }
    return v;
}

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
    B(0,N-1) = beta;
    B(N-1,0) = beta;

    // Build the diagonal matrix which is on the subdiagonals
    MatrixXd C = MatrixXd::Zero(N,N);
    C.diagonal() = beta*VectorXd::Ones(N);

    // Populate A with the blocks
    MatrixXd A = MatrixXd::Zero(N*N,N*N);
    for(int i=0; i<N; i++)
    {
        A.block(i*N,i*N,N,N) = B;
        if(i==0)
            A.block((N-1)*N,i*N,N,N) = C;
        else
            A.block((i-1)*N,i*N,N,N) = C;
        if(i==N-1)
            A.block(0,i*N,N,N) = C;
        else
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
    B.diagonal(1) = VectorXd::Ones(N-1);
    B.diagonal(-1) = VectorXd::Ones(N-1);
    B(0,N-1) = 1;
    B(N-1,0) = 1;

    // Build the diagonal matrix which is on the subdiagonals
    MatrixXd C = MatrixXd::Identity(N,N);

    // Populate A with the blocks
    MatrixXd A = MatrixXd::Zero(N*N,N*N);
    for(int i=0; i<N; i++)
    {
        A.block(i*N,i*N,N,N) = B;
        if(i==0)
            A.block((N-1)*N,i*N,N,N) = C;
        else 
            A.block((i-1)*N,i*N,N,N) = C;
        if(i==N-1)
            A.block(0,i*N,N,N) = C;
        else
            A.block((i+1)*N,i*N,N,N) = C;
    }

    return A;
}

void updateLoadU(VectorXd u, VectorXd v, int N, double dt, 
        double Re, VectorXd &f_U)
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
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            // Point to be evaluated 
            u_0 = u[i*N+j];

            // x-velocity to the north
            if(j==N-1) u_N = u[i*N]; // north boundary
            else u_N = u[i*N+(j+1)];

            // x-velocity to the east
            if(i==N-1) u_E = u[j]; // east boundary
            else u_E = u[(i+1)*N+j];

            // x-velocity to the south
            if(j==0) u_S = u[i*N+N-1]; // south boundary
            else u_S = u[i*N+(j-1)];

            // x-velocity to the west
            if(i==0) u_W = u[(N-1)*N+j]; // west boundary
            else u_W = u[(i-1)*N+j];

            // y-velocity north-west
            if(i==0)
            {
                if(j==N-1) v_NW = v[N-1]; // NW corner
                else v_NW = v[(j+1)*N+(N-1)]; // west border
            }
            else if(j==N-1) v_NW = v[i-1]; // north border
            else v_NW = v[(j+1)*N+(i-1)];

            // y-velocity north-east
            if(j==N-1) v_NE = v[i]; // north border
            else v_NE = v[(j+1)*N+i];

            // y-velocity south-east
            v_SE = v[j*N+i];

            // y-velocity south-west
            if(i==0) v_SW = v[j*N+(N-1)]; // west border
            else v_SW = v[j*N+(i-1)];
          
            // compute load vector
            f_U[i*N+j] = u_0-alpha*(
                    u_0*(u_E-u_W)+0.5*u_0*(v_NW+v_NE-v_SW-v_SE)
                    +0.125*(u_N-u_S)*(v_NW+v_NE+v_SE+v_SW));
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
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            // Point to be evaluated 
            v_0 = v[i*N+j];

            // y-velocity to the east
            if(j==N-1) v_E = v[i*N]; // east boundary
            else v_E = v[i*N+(j+1)];

            // y-velocity to the north
            if(i==N-1) v_N = v[j]; // north boundary
            else v_N = v[(i+1)*N+j];

            // y-velocity to the west
            if(j==0) v_W = v[i*N+N-1]; // west boundary
            else v_W = v[i*N+(j-1)];

            // y-velocity to the south
            if(i==0) v_S = v[(N-1)*N+j]; // west boundary
            else v_S = v[(i-1)*N+j];

            // x-velocity south-east
            if(i==0)
            {
                if(j==N-1) u_SE = u[N-1]; // SE corner
                else u_SE = u[(j+1)*N+(N-1)]; // south border
            }
            else if(j==N-1) u_SE = u[i-1]; // east border
            else u_SE = u[(j+1)*N+(i-1)];

            // x-velocity north-east
            if(j==N-1) u_NE = u[i]; // east border
            else u_NE = u[(j+1)*N+i];

            // x-velocity north-west
            u_NW = u[j*N+i];

            // x-velocity south-west
            if(i==0) u_SW = u[j*N+(N-1)]; // south border
            else u_SW = u[j*N+(i-1)];

            // compute load vector
            f_V[i*N+j] = v_0-alpha*(
                    v_0*(v_N-v_S)+0.5*v_0*(u_NE+u_SE-u_NW-u_SW)
                    +0.125*(v_E-v_W)*(u_NW+u_NE+u_SE+u_SW));
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
            // x-velocity west of pressure point
            U_W = U[i*N+j];

            // x-velocity east, periodic on eastern border
            if(i==N-1) U_E = U[j];
            else U_E = U[(i+1)*N+j];

            // y-velocity south
            V_S = V[j*N+i];

            // y-velocity north, periodic on northern border
            if(j==N-1) V_N = V[i];
            else V_N = V[(j+1)*N+i];
          
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
    double p_E, p_W, p_N, p_S; 

    // Approximation to the x-derivative of pressure
    VectorXd pDiffX(N*N);

    // Approximation to the y-deriavtive of pressure
    VectorXd pDiffY(N*N);

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            p_E = p[i * N + j];
            if( i == 0 ) 
                p_W = p[( N - 1 ) * N + j];
            else 
                p_W = p[( i - 1 ) * N + j];

            p_N = p[j * N + i];
            if( i == 0 ) 
                p_S = p[j * N + ( N - 1 )];
            else 
                p_S = p[j * N + ( i - 1 )];

            pDiffX[i*N+j] = p_E - p_W; 
            pDiffY[i*N+j] = p_N - p_S; 
        }
    }

    // Compute
    u = U - dt*N*pDiffX;
    v = V - dt*N*pDiffY;
}

VectorXd buildVorticityVector(VectorXd u, VectorXd v, int N)
{
    double u_N, u_S, v_W, v_E; 
    VectorXd omega(N*N);
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            u_N = u[i*N+j];

            if(j==0) u_S = u[i*N+(N-1)];
            else u_S = u[i*N+(j-1)];

            v_E = v[j*N+i];

            if(i==0) v_W = v[j*N+(N-1)];
            else v_W = v[j*N+(i-1)];

            omega[i*N+j] = u_N-u_S-(v_E-v_W);
        }
    }
    return omega*N;
}

#endif
