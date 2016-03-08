#include "Functions.h"

int main(int argc, char* argv[]) 
{
    int N; 
    double Re; 

    if(argc==1)
    {
        cout << "N: ";
        cin >> N; 
        cout << "Re: ";
        cin >> Re; 
    }
    else
    {
        N = atoi(argv[1]);      // Number of cells per direction
        Re = atof(argv[2]);     // Reynolds number
    }
    if(N==20) N=21;             //TODO: FInd out why N=20 doesnt

    int i, j;                   // Counters for loops
    double a = 1.0;             // Velocity of lid
    double dx = 1.0/N;          // Grid spacing
    double t = 0.0;             // Time counter
    double dt = 500*dx/Re;      // Time step
    if(Re>1000) dt /= 5.0;      // Ensure time step is small enough for high-Re
    double tol = 1e-5;          // Tolerance for convergence

    ofstream fs;           // File stream for writing res
    char filename[50];

    // Initiate solution vectors
    VectorXd u = VectorXd::Zero(N*(N-1));
    VectorXd prev = u;
    VectorXd eps = VectorXd::Ones(N*(N-1));
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

    cout << "iter" << "\t" << "time" << "\t" << "||u_c||" << "\t" << "relCha" << "\n";

    int iter = 0;
    while(eps.lpNorm<Infinity>()/u.lpNorm<Infinity>()>tol)
    {
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

        eps = u-prev;
        prev = u;

//        if(iter%10==0)
            cout << iter << "\t" << iter*dt << "\t" 
                << eps.squaredNorm() << "\t" << eps.lpNorm<Infinity>()/u.lpNorm<Infinity>() << "\n";
    }

    cout << "Steady state reached." << endl;

    // Calculate vorticity vector and stream function
    VectorXd omega = buildVorticityVector(u, v, N);
    VectorXd psi = solver_psi.solve(omega);

    // Print results to file

    sprintf(filename, "../Results/psi_N%d_Re%d.out", N, int(Re));
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

    sprintf(filename, "../Results/omega_N%d_Re%d.out", N, int(Re));
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

    sprintf(filename, "../Results/uGC_N%d_Re%d.out", N, int(Re));
    fs.open(filename, std::fstream::out);
    i= N/2-1;
    for(j=0; j<N; j++)
    {
        fs << (j+0.5)*dx << "\t" << u[i*N+j] << "\n";
    }
    fs.close();

    sprintf(filename, "../Results/vGC_N%d_Re%d.out", N, int(Re));
    fs.open(filename, std::fstream::out);
    for(j=0; j<N; j++)
    {
        fs << (j+0.5)*dx << "\t" << v[i*N+j] << "\n";
    }
    fs.close();

    sprintf(filename, "../Results/p_N%d_Re%d.out", N, int(Re));
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

    sprintf(filename, "../Results/vec_N%d_Re%d.out", N, int(Re));
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
}
