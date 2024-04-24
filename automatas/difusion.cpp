#include <iostream>
#include <cmath>
#include "Random64.h"

const int Lx=1024;
const double p=0.5;

const int Q=2;

class lattice_gas
{
private:
    //V[0]:derecha , V[1]:izquierda 
    int V[Q];
    int n[Lx][Q];
    int nnew[Lx][Q];
public:
    lattice_gas(void);
    void borrese(void);
    void inicie(int N, double mu, double sigma, Crandom & ran64);
    void show(void);
    void colisione(Crandom & ran64);
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
    {
        for (int i = 0; i < Q; i++)
        {
            n[ix][i]=0;
        }
    }
    
}

void lattice_gas::inicie(int N, double mu, double sigma, Crandom & ran64)
{
    int ix, i;
    while (N>0)
    {
        //Escoge al azar un lugar y una direcion
        ix = (int) ran64.gauss(mu, sigma);
        if(ix<0) ix=0;
        if(ix>=Lx) ix=Lx-1;

        if (ran64.r()<0.5) i=0; else i=1;       

        if(n[ix][i]==0)
        {
            n[ix][i]++;
            N--;
        }
    }
}

void lattice_gas::show(void)
{
    for (int i = 0; i < Q; i++)
    {
        for (int ix = 0; ix < Lx; ix++)
        {
            std::cout<<n[ix][i]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void lattice_gas::colisione(Crandom & ran64)
{
    for (int ix = 0; ix < Lx; ix++)
    {
        if (ran64.r()<p) 
            for (int i = 0; i < Q; i++)
                nnew[ix][i]=n[ix][i]; //lo dejo igual

        else
            for (int i = 0; i < Q; i++)
                nnew[ix][i]=n[ix][(i+1)%2]; //intercambio los contenidos

    }
    
}

void lattice_gas::adveccione(void)
{
    for ( int ix = 0; ix < Lx; ix++)
        for (int i = 0; i < Q; i++)
        {
            n[(ix+V[i]+Lx)%Lx][i] = nnew[ix][i]; //se copia el valor nuevo
        }
}

double lattice_gas::rho(int ix, bool UseNew)
{
    if(UseNew)
        return nnew[ix][0] + nnew[ix][1];
    else    
        return n[ix][0] + n[ix][1];
}

double lattice_gas::get_sigma2(void)
{
    double prom=0, sigma2=0;
    int ix, N=0;
    //Calcular el promedio
    for (int ix = 0; ix < Lx; ix++)
    {
        N += rho(ix, false);
        prom += ix*rho(ix, false);
    }
    prom/=N;

    //Calcular sigma2
    for(ix=0;ix<Lx;ix++)
        sigma2+=pow(ix-prom,2.0)*rho(ix,false);
    sigma2/=(N-1);

    return sigma2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main (void)
{
    lattice_gas difusion;
    Crandom ran64(1);

    int N = (int) Lx*0.1;
    double mu= Lx/2;
    double sigma= Lx/8;
    int t, tmax=400;

    difusion.borrese();
    difusion.inicie(N, mu, sigma, ran64);

    for ( t = 0; t < tmax; t++)
    {
        std::cout<<t<<" "<<difusion.get_sigma2()<<std::endl;
        difusion.colisione(ran64);
        difusion.adveccione();
    }

    return 0;
}