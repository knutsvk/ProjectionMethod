#include "Functions.h"

int main(void) 
{
    int N = 5;             // Number of cells per direction
    
    double Re = 100.0;      // Reynolds number
    double a = 1.0;         // Velocity of lid
    double tau = 0.1/N;     // End time TODO: run until steady?
    double dx = 1.0/N;      // Grid spacing
    double t = 0.0;         // Time counter
    double dt = dx/10.0;    // Time step TODO: check stability

    // Initiate solution vectors
    VectorXd u = VectorXd::Zero(N*(N-1));
    VectorXd v = VectorXd::Zero(N*(N-1));
    VectorXd p = VectorXd::Zero(N*N);

    // Initiate intermediate velocity vectors
    VectorXd U = VectorXd::Zero(N*(N-1));
    VectorXd V = VectorXd::Zero(N*(N-1));

    // Declare load vectors 
    VectorXd f_U(N*(N-1));
    VectorXd f_V(N*(N-1));
    VectorXd f_p(N*N);

    // TODO: Make the matrices sparse
    // Build matrix for velocity equations
    MatrixXd A_UV = buildVelocityMatrix(N, dt, Re);
    cout << "A_UV = " << endl << A_UV << endl; 

    // Build matrix for pressure equations
    MatrixXd A_p = buildPressureMatrix(N); 
    cout << "A_p = " << endl << A_p << endl; 

    while(t<tau)
    {
        // Make sure we don't go past end time
        if(t+dt>tau) dt = tau-t;

        // Compute rhs for update formula for x-velocity
        updateLoadU(u, v, N, dt, a, Re, f_U);
        cout << "f_U = " << endl << f_U << endl; 

        // Solve Au * U = bu
        
        // Repeat the above for y-velocity
        updateLoadV(u, v, N, dt, a, Re, f_V);
        cout << "f_V = " << endl << f_V << endl; 
        
        // Compute rhs for update formula for pressure
        // Solve Ap * p = bp
        
        // Find u and v at next time step
        
        t += dt;
    }

    // Print results (remember ghost points!) to file
}
