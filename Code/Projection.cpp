#include "Functions.h"

int main(int argc, char* argv[]) 
{
    int N = atoi(argv[1]);     // Number of cells per direction
    int i, j;
  
    double Re = atof(argv[2]); // Reynolds number
    double a = 1.0;            // Velocity of lid
    double dx = 1.0/N;         // Grid spacing
    double t = 0.0;            // Time counter
    double dt = dx/1.0;       // Time step TODO: check stability
    double tau = 1000*dt;      // End time TODO: run until steady?

    ofstream fs;           // File stream for writing res
    char filename[20];

    // Initiate solution vectors
    VectorXd u = VectorXd::Zero(N*(N-1));
    VectorXd v = VectorXd::Zero(N*(N-1));
    VectorXd p = VectorXd::Ones(N*N);

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
        updateLoadV(u, v, N, dt, Re, f_V);
        V = solver_UV.solve(f_V);
        
        // Compute rhs for update formula for pressure
        updateLoadp(U, V, N, dt, f_p);

        // Solve Ap * p = bp
        p = solver_p.solve(f_p);
        
        // Find u and v at next time step
        updateVelocities(U, V, p, N, dt, u, v);
        
        t += dt;
        iter++; 
        if(iter%10==0)
        {
            cout << "iteration: " << iter << endl; 
            sprintf(filename, "../Results/uGC%d.out",iter/10);
            fs.open(filename);
            i= N/2-1;
            for(j=0; j<N; j++)
            {
                fs << (i+1)*dx << "\t" << (j+0.5)*dx << "\t" 
                    << u[i*N+j] << "\n";
            }
            fs.close();
        }
    }

    cout << "Simulation complete!" << endl; 

    // Print results (TODO: remember ghost points?) to file
    fs.open("../Results/uGC.out", std::fstream::out);
    i= N/2-1;
    for(j=0; j<N; j++)
    {
        fs << (i+1)*dx << "\t" << (j+0.5)*dx << "\t" 
            << u[i*N+j] << "\n";
    }
    fs.close();

    fs.open("../Results/vGC.out", std::fstream::out);
    for(j=0; j<N; j++)
    {
        fs << (j+0.5)*dx << "\t" << (i+1)*dx << "\t" 
            << v[i*N+j] << "\n";
    }
    fs.close();

    fs.open("../Results/p.out", std::fstream::out);
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            fs << (i+1)*dx << "\t" << (j+1)*dx << "\t" 
                << p[i*N+j] << "\n";
        }
        fs << "\n";
    }
    fs.close();

    fs.open("../Results/uABS.out", std::fstream::out);
    for(i=0; i<N-1; i++)
    {
        for(j=0; j<N-1; j++)
        {
            fs << (i+1)*dx << "\t" << (j+1)*dx << "\t" 
                << sqrt(pow(u[i*N+j],2)+pow(v[j*N+i],2)) << "\n";
        }
        fs << "\n";
    }
    fs.close();

    fs.open("../Results/vec.out", std::fstream::out);
    for(i=0; i<N-1; i++)
    {
        for(j=0; j<N-1; j++)
        {
            fs << (i+1)*dx << "\t" << (j+1)*dx << "\t" 
                << u[i*N+j] << "\t" << v[j*N+i] << "\n";
        }
    }
    fs.close();
}
