#include <iostream>
#include <math.h>

using namespace std;

int factorial(int x)
{
    //Devulve el factorial de x
    if (x==0)
        return 1;
    else
        return x*factorial(x-1);
}

int main()
{
    int x = 5;
    cout<<factorial(x)<<endl;
    return 0;
}
