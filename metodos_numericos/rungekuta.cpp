#include <iostream>
#include <math.h>

using namespace std;

double f(double t, double x)
{
    return x;
}

void unpasoRK(double & t, double & x, double dt)
{
    double dx1, dx2, dx3, dx4;
    dx1 = dt*f(t,x);
    dx2 = dt*f(t+dt/2,x+dx1/2);
    dx3 = dt*f(t+dt/2,x+dx2/2);
    dx4 = dt*f(t+dt,x+dx3);    
    t+=dt; x+=(dx1+2*(dx2+dx3)+dx4)/6;
}

int main()
{
    double t, x, dt=0.0625;

    for (t=0,x=1; t <= 2; )
    {
        cout<<t<<" "<<x<<" "<<exp(t)<<endl;
        unpasoRK(t,x,dt);
    }
    return 0;
}