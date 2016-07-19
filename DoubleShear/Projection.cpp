#include "Functions.h"

int main(int argc, char* argv[]) 
{
    int N; 

    if(argc==1)
    {
        cout << "N: ";
        cin >> N; 
    }
    else
    {
        N = atoi(argv[1]);      // Number of cells per direction
    }
    if(N==20) N=21;      

    int i, j;                   // Counters for loops
    double Re = 30.0 / 2.0e-4;  // Reynolds number
    double dx = 1.0/N;          // Grid spacing
    double t = 0.0;             // Time counter
    double dt = 500*dx/Re;      // Time step
    double tStop = 2.0;
    if(Re>1000) dt /= 5.0;      // Ensure clock step is small enough for high-Re

    ofstream fs;           // File stream for writing res
    char filename[50];

    cout << "Incompressible flow in a lid-driven cavity." << endl
        << "Runtime parameters: N = " << N << ", Re = " << Re 
        << endl << "Building matrices... " << endl; 

    clock_t initiationStart = clock(); 

    // Initiate solution vectors
    VectorXd u = VectorXd::Zero(N*(N-1));
    VectorXd v = VectorXd::Zero(N*(N-1));
    VectorXd p(N*N);

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
    // Get factorization
    SimplicialLDLT<SparseMatrix<double> > solver_UV;
    solver_UV.compute(Sp_A_UV);

    // Build matrix for pressure equations, get solver
    MatrixXd A_p = buildPressureMatrix(N); 
    SparseMatrix<double> Sp_A_p = A_p.sparseView();
    SimplicialLDLT<SparseMatrix<double> > solver_p;
    solver_p.compute(Sp_A_p);

    // Build matrix for pressure equations, get solver
    MatrixXd A_psi = buildStreamMatrix(N-1); 
    SparseMatrix<double> Sp_A_psi = A_psi.sparseView();
    SimplicialLDLT<SparseMatrix<double> > solver_psi;
    solver_psi.compute(Sp_A_psi);

    clock_t initiationEnd = clock(); 

    cout << "Building of matrices complete. Time elapsed: " 
        << double(initiationEnd-initiationStart) / CLOCKS_PER_SEC 
        << " seconds." << endl 
        << "Advancing in time until final time..." << endl; 

    clock_t iterateStart = clock(); 

    int iter = 0;
    while(t < tStop)
    {
        // Make sure we dont go past tStop
        if( t + dt > tStop )
            dt = tStop - t; 

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
        
        // Find u and v at next clock step
        updateVelocities(U, V, p, N, dt, u, v);
        
        t += dt;
        iter++; 
    }

    clock_t iterateEnd = clock();

    cout << "Final time reached after " << iter 
        << " iterations. " << "Time elapsed: " 
        << double(iterateEnd-iterateStart) / CLOCKS_PER_SEC 
        << " seconds." << endl << "Time per iteration: " 
        << double(iterateEnd-iterateStart) / CLOCKS_PER_SEC / iter
        << " seconds." << endl; 

    // Calculate vorticity vector and stream function
    VectorXd omega = buildVorticityVector(u, v, N);
    VectorXd psi = solver_psi.solve(omega);

    cout << "Printing results to file... " << endl; 

    // Print results to file

    sprintf(filename, "../Results/psi_N%d.out", N);
    fs.open(filename, std::fstream::out);
    for(i=0; i<N-1; i++)
    {
        for(j=0; j<N-1; j++)
        {
            fs << (i+1)*dx << "\t" << (j+1)*dx << "\t" 
                << psi[i*(N-1)+j] << "\n";
        }
        fs << "\n";
    }
    fs.close();

    sprintf(filename, "../Results/omega_N%d.out", N);
    fs.open(filename, std::fstream::out);
    for(i=0; i<N-1; i++)
    {
        for(j=0; j<N-1; j++)
        {
            fs << (i+1)*dx << "\t" << (j+1)*dx << "\t" 
                << omega[i*(N-1)+j] << "\n";
        }
        fs << "\n";
    }
    fs.close();

    sprintf(filename, "../Results/uGC_N%d.out", N);
    fs.open(filename, std::fstream::out);
    i= N/2-1;
    for(j=0; j<N; j++)
    {
        fs << (j+0.5)*dx << "\t" << u[i*N+j] << "\n";
    }
    fs.close();

    sprintf(filename, "../Results/vGC_N%d.out", N);
    fs.open(filename, std::fstream::out);
    for(j=0; j<N; j++)
    {
        fs << (j+0.5)*dx << "\t" << v[i*N+j] << "\n";
    }
    fs.close();

    sprintf(filename, "../Results/p_N%d.out", N);
    fs.open(filename, std::fstream::out);
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

    sprintf(filename, "../Results/vec_N%d.out", N);
    fs.open(filename, std::fstream::out);
    fs << "x\ty\tu\tv\n";
    for(i=0; i<N-1; i++)
    {
        for(j=0; j<N-1; j++)
        {
            fs << (i+1)*dx << "\t" << (j+1)*dx << "\t" 
                << u[i*N+j] << "\t" << v[j*N+i] << "\n";
        }
    }
    fs.close(); 

    cout << "All done. Total runtime: " 
        << double(clock() - initiationStart) / CLOCKS_PER_SEC 
        << " seconds. " << endl; 
}
