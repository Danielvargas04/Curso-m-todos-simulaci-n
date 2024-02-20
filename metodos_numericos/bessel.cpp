#include <iostream>
#include <math.h>

using namespace std;

double f(double alpha, double x, double t)
{
    return cos(alpha*t-x*sin(t));
}

double int_simpson(double alpha, double x, double a, double b, int n )
{
    double h, t, suma=0;
    n*=2;
    h = (b-a)/n;
    

    for(int i = 0; i<=n; i++)
    {
        t=a + i*h;
        if(i==0 || i==n)
            suma+=f(alpha, x, t);
        else if (i%2==0)
            suma+=2*f(alpha, x, t);
        else
            suma+=4*f(alpha, x, t);
        
    }
    return suma*h/3;
}

double bessel(double alpha, double x)
{
    double a=0, b=M_PI;
    int n=50;
    return int_simpson(alpha, x, a, b, n)/M_PI;
}

int main()
{
    /*Funcion de bessel para alpha = 3*/
    double alpha=3, x=0.1;
    for (double  x = 0; x <=10; x+=0.1)
    {
        cout<<x<<" "<<bessel(alpha, x)<<endl;
    }
    
    return 0;
}