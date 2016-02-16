#include "Functions.h"

int main(void) 
{
    int N =50;             // Number of cells per direction
    
    double Re = 100.0;      // Reynolds number
    double a = 1.0;         // Velocity of lid
    double tau = 1.0;     // End time TODO: run until steady?
    double dx = 1.0/N;      // Grid spacing
    double t = 0.0;         // Time counter
    double dt = dx/5.0;    // Time step TODO: check stability

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

    // TODO: Make the matrices sparse from beginning
    // Build matrix for velocity equations
    MatrixXd A_UV = buildVelocityMatrix(N, dt, Re);
    SparseMatrix<double> Sp_A_UV = A_UV.sparseView();
    cout << "Finished building matrix A_UV!" << endl; 

    // Get factorization
    SimplicialLDLT<SparseMatrix<double> > solver_UV;
    solver_UV.compute(Sp_A_UV);
    if(solver_UV.info()!=Success)
        cout << "Decomposition of matrix A_UV failed!" << endl; 
    else
        cout << "Finished decomposition of matrix A_UV" << endl; 

    // Build matrix for pressure equations, get solver
    MatrixXd A_p = buildPressureMatrix(N); 
    SparseMatrix<double> Sp_A_p = A_p.sparseView();
    cout << "Finished building matrix A_p!" << endl; 
    SimplicialLDLT<SparseMatrix<double> > solver_p;
    solver_p.compute(Sp_A_p);
    if(solver_p.info()!=Success)
        cout << "Decomposition of matrix A_p failed!" << endl; 
    else
        cout << "Finished decomposition of matrix A_p" << endl; 

    int iter = 0;
    while(t<tau)
    {
        // Make sure we don't go past end time
        if(t+dt>tau) dt = tau-t;

        // Compute rhs for update formula for x-velocity
        updateLoadU(u, v, N, dt, a, Re, f_U);

        // Solve Au * U = bu
        U = solver_UV.solve(f_U);
        
        // Repeat the above for y-velocity
        updateLoadV(u, v, N, dt, a, Re, f_V);
        V = solver_UV.solve(f_V);
        
        // Compute rhs for update formula for pressure
        updateLoadp(u, v, N, dt, f_p);

        // Solve Ap * p = bp
        p = solver_p.solve(f_p);
        
        // Find u and v at next time step
        updateVelocities(U, V, p, N, dt, u, v);
        
        t += dt;
        iter++; 
        if(iter%100==0) cout << "iteration: " << iter << endl; 
    }

    cout << "Simulation complete!" << endl; 

    // Print results (TODO: remember ghost points?) to file
    std::fstream fs; 
    fs.open("../Results/uGC.out", std::fstream::out);
    int i = N/2-1;
    for(int j=0; j<N; j++)
    {
        fs << (i+1)*dx << "\t" << (j+0.5)*dx << "\t" 
            << u[i*N+j] << "\n";
    }
    fs.close();

    fs.open("../Results/vGC.out", std::fstream::out);
    for(int j=0; j<N; j++)
    {
        fs << (j+0.5)*dx << "\t" << (i+1)*dx << "\t" 
            << v[i*N+j] << "\n";
    }
    fs.close();

}
