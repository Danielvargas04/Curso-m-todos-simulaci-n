#include <iostream>
#include <cmath>
#include "Random64.h"

const int Lx=128;
const double p=0.5;

const int Q=2;

class lattice_gas
{
private:
    //V[0]:derecha , V[1]:izquierda 
    int V[Q];
    double f[Lx][Q],fnew[Lx][Q]; // f[ix][i]
public:
    lattice_gas(void);
    void borrese(void);
    void inicie(double mu, double sigma);
    void show(void);
    void colisione(void);
    void adveccione(void);
    double rho(int ix, bool UseNew);
    double get_sigma2(void);
};

//--------- Metodos de la clase lattice_gas---------

lattice_gas::lattice_gas(void)
{
    //definir los vectores velocidad
    V[0]=1; V[1]=-1;
}

void lattice_gas::borrese(void)
{
    for (int ix = 0; ix < Lx; ix++)
        for (int i = 0; i < Q; i++)
        {
            f[ix][i]=0;
        }    
}

void lattice_gas::inicie(double mu, double sigma)
{
    for(int ix=0;ix<Lx;ix++)
        for(int i=0;i<Q;i++)
            f[ix][i]=fnew[ix][i]=0.5/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2.0));
}

void lattice_gas::show(void)
{
    for(int ix=0;ix<Lx;ix++)
        std::cout<<ix<<" "<<rho(ix,true)<<std::endl;
    std::cout<<std::endl;
}

void lattice_gas::colisione(void)
{
    int ix,i,j;
    for(ix=0;ix<Lx;ix++) //para cada celda
        for(int i=0;i<Q;i++)
        {
            j=(i+1)%2; 
            fnew[ix][i]=p*f[ix][i]+(1-p)*f[ix][j];
        }    
}

void lattice_gas::adveccione(void)
{
    for ( int ix = 0; ix < Lx; ix++)
        for (int i = 0; i < Q; i++)
        {
            f[(ix+V[i]+Lx)%Lx][i] = fnew[ix][i]; //se copia el valor nuevo
        }
}

double lattice_gas::rho(int ix, bool UseNew)
{
    if(UseNew)
        return fnew[ix][0] + fnew[ix][1];
    else    
        return f[ix][0] + f[ix][1];
}

double lattice_gas::get_sigma2(void)
{
    double prom=0, sigma2=0, N=0;
    //Calcular el promedio
    for (int ix = 0; ix < Lx; ix++)
    {
        N += rho(ix, true);
        prom += ix*rho(ix, true);
    }
    prom/=N;

    //Calcular sigma2
    for(int ix=0; ix < Lx; ix++)
        sigma2+=pow(ix-prom,2.0)*rho(ix,true);
    sigma2/=N;

    return sigma2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main (void)
{
    lattice_gas difusion;
    Crandom ran64(1);

    int N = 1;
    double mu= Lx/2;
    double sigma= Lx/64;
    int t, tmax=400;

    difusion.borrese();
    difusion.inicie(mu, sigma );

    for ( t = 0; t < tmax; t++)
    {
        std::cout<<t<<" "<<difusion.get_sigma2()<<std::endl;
        difusion.colisione();
        difusion.adveccione();
    }
    //difusion.show();
    
    return 0;
}