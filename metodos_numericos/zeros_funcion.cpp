#include <iostream>
#include <math.h>

using namespace std;

double f(double x){
    return sin(x)/x;
}

const double Tol=1e-7;

int main(){

    double a=3, b=4, fa=f(a), fm;
    double m;

    while (b-a > Tol){

        m = (b+a)/2;
        fm = f(m);
        if (fa*fm > 0)
        {
            a = m;
            fa = fm;
        }
        else
            b = m;
    }

    cout<<"El cero es "<<(b+a)/2<<endl;

    return 0;
}