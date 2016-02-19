#include "Functions.h"

int main(void)
{
    int N=5;
    VectorXd u(N*N);
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            u[i*N+j] = i*N+(j+1);
        }
    }
    cout << "N=" << N << endl;
    cout << "u=" << endl << u << endl;
    cout << "u.tail(N*(N-1))=" << endl << u.tail(N*(N-1)) << endl;
    cout << "u.head(N*(N-1))=" << endl << u.head(N*(N-1)) << endl;
    MatrixXd A = buildPressureMatrix(N);
    cout << "A=" << endl << A << endl;
}
