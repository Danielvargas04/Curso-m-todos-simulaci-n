#include <iostream>
#include <math.h>

using namespace std;

double f(double t, double x)
{
    return x;
}

void unpasodeuler(double &t, double &x, double dt)
{
    double dx;
    dx = dt*f(t,x);
    t+=dt; x+=dx;
}

int main()
{
    double t, x, dt=0.01;

    for (t=0,x=1; t <= 2; )
    {
        cout<<t<<" "<<x<<" "<<exp(t)<<endl;
        unpasodeuler(t,x,dt);
    }
    return 0;
}