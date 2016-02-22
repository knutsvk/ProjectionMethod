#include "Functions.h"

int main(void)
{
    double dx=.01;
    double x=0.7349;
    int i=round(x/dx)-1;
    x=(i+1)*dx;
    cout << x << endl;

    double u, v;
    for(int i=0; i<11; i++)
    {
        u=1-2.0*i/10.0;
        for(int j=0; j<11; j++)
        {
            v=1-2.0*j/10.0;
            cout << "u= " << u << "\t" 
                << "v= " << v << "\t"
                << "angle= " << atan(v/u)+M_PI/2.0 << endl; 
        }
    }
}
