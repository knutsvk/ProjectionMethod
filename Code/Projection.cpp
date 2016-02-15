#include "Functions.h"

int main(void) 
{
    int N = 5;             // Number of cells per direction
    
    double Re = 100.0;      // Reynolds number
    double a = 1.0;         // Velocity of lid
    double tau = 1.0;       // End time
    double dx = 1.0/N;      // Grid spacing
    double t = 0.0;         // Time counter
    double dt = dx/10.0;    // Time step

    // Initiate solution vectors
    VectorXd u = VectorXd::Zero(N*(N-1));
    VectorXd v = VectorXd::Zero(N*(N-1));
    VectorXd p = VectorXd::Zero(N*N);

    // Initiate intermediate velocity vectors
    VectorXd U = VectorXd::Zero(N*(N-1));
    VectorXd V = VectorXd::Zero(N*(N-1));

    // Build matrix for velocity equations
    MatrixXd Au = buildVelocityMatrix(N, dt, Re);
    cout << "Au = " << endl << Au << endl; 

    // Build matrix for pressure equations
    MatrixXd Ap = buildPressureMatrix(N); 
    cout << "Ap = " << endl << Ap << endl; 
    

    while(t<tau)
    {
        // Make sure we don't go past end time
        if(t+dt>tau) dt = tau-t;

        // Compute rhs for update formula for x-velocity
        // Add BCs to rhs 
        // Solve Au * U = bu
        
        // Repeat the above for y-velocity
        
        // Compute rhs for update formula for pressure
        // Solve Ap * p = bp
        
        // Find u and v at next time step
        
        t += dt;
    }
}
