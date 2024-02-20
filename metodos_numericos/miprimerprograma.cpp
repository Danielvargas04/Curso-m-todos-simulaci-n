#include <iostream>
#include <math.h>

using namespace std;

int main(){

    for (double x = 0.1; x <= 10; x+=0.1)
    {
        cout<<x<<" "<<sin(x)/x<<endl;
    }
    
    return 0;
}