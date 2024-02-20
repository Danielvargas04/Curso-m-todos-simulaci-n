#include <iostream>
#include <math.h>

using namespace std;

double f(double x)
{
    return cos(-0.1*sin(x));
}

int main()
{
    double x, a=0, b=M_PI, h, suma=0;
    int n = 100;
    h = (b-a)/n;

    for(int i = 0; i<=n; i++)
    {
        x=a + i*h;
        if(i==0 || i==n)
            suma+=f(x);
        else if (i%2==0)
            suma+=2*f(x);
        else
            suma+=4*f(x);
        
    }
    cout<<"La intergal es "<< suma*h/(3*M_PI)<<endl;
    return 0;
}