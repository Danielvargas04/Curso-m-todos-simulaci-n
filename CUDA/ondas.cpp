#include <iostream>
#include <fstream>
#include <cmath>

//Constantes del problemas

const int Lx=128;
const int Ly=128;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; //C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//---------------class Laticce--------------
class LatticeBoltzmann
{
private:
    double w[Q];            //Weights
    int Vx[Q], Vy[Q];    //Velocity vectors
    double *f, *fnew;            //Weigtha
public:
    LatticeBoltzmann(void);
    ~LatticeBoltzmann();
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    //------------ Campos macroscopicos-----------------
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    //------------ Funciones de equilibrio-----------------
    double feq(double rho0, double Jx0, double Jy0, int i);
    //---------- Evolucion temporal----------------
    void Star(double rho0, double Jx0, double Jy0);
    void Collision(void);
    void ImposeFields(int t);
    void Adveccion(void);
    //---------- Funciones Globales----------------
    void Print(const  char * NameFile);
};

//-------------Implementacion de funciones-----------------

LatticeBoltzmann::LatticeBoltzmann(void)
{
    //Set the weights
    w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
    //Set Velocity vectors
    Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
    Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;
    //Create the dynamic arrays
    int ArraySize=Lx*Ly*Q;
    f=new double [ArraySize]; fnew=new double [ArraySize];
}

LatticeBoltzmann::~LatticeBoltzmann()
{
    delete[] f; delete[] fnew;
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew)
{
    double sum=0.0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        if (UseNew) sum+=fnew[n0];
        else sum+=f[n0];
    }
    return sum;    
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew)
{
    double sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        if (UseNew) sum+=Vx[i]*fnew[n0];
        else sum+=Vx[i]*f[n0];
    }
    return sum;    
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew)
{
    double sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        if (UseNew) sum+=Vy[i]*fnew[n0];
        else sum+=Vy[i]*f[n0];
    }
    return sum;    
}

double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i)
{
    if(i>0)
        return 3*w[i]*(C2*rho0 + Vx[i]*Jx0 + Vy[i]*Jy0);
    else
        return rho0*AUX0;
}

void LatticeBoltzmann::Star(double rho0, double Jx0, double Jy0)
{
    int ix, iy, i, n0;
    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
            for(i = 0; i < Q; i++)  //on each direction
            {
                n0=n(ix,iy,i);
                f[n0]=feq(rho0, Jx0, Jy0, i);
            }
}

void LatticeBoltzmann::Collision(void)
{
    int ix, iy, i, n0;
    double rho0, Jx0, Jy0;

    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
        {
            rho0 = rho(ix, iy, false);
            Jx0 = Jx(ix, iy, false);
            Jy0 = Jy(ix, iy, false);
            for(i = 0; i < Q; i++)  //on each direction
            {
                n0=n(ix,iy,i);
                fnew[n0]=UmUtau*f[n0] + Utau*feq(rho0, Jx0, Jy0, i);
            }
        }
}

void LatticeBoltzmann::ImposeFields(int t)
{
    int ix, iy, i, n0;
    double rho0, Jx0, Jy0;
    double lambda=10, omega=2*M_PI/lambda*C;

    //An oscillating source in the middle
    ix=Lx/2; iy=Ly/2;
    rho0=10*sin(omega*t); Jx0=Jx(ix, iy, false); Jy0=Jy(ix, iy, false);
    for (i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        fnew[n0] = feq(rho0, Jx0, Jy0, i);
    }
}

void LatticeBoltzmann::Adveccion(void)
{
    int ix, iy, i, ixnew, iynew, n0, n0new;
    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
            for(i = 0; i < Q; i++)  //on each direction
            {
                ixnew=(ix+Vx[i]+Lx)%Lx; iynew=(iy+Vy[i]+Ly)%Ly;
                n0=n(ix, iy, i); n0new=n(ixnew, iynew, i);
                f[n0new]=fnew[n0];  //periodic boundaries
            }
}



//--------------Funciones Globales--------------

void LatticeBoltzmann::Print(const char * NameFile)
{
    std::ofstream MyFile(NameFile); double rho0; int ix,iy;
    for (ix = 0; ix < Lx; ix++)
    {
        for (iy = 0; iy < Ly; iy++)
        {
            rho0=rho(ix, iy, true);
            MyFile<<ix<<" "<<iy<<" "<<rho0<<std::endl;
        }
        MyFile<<std::endl;
    }
    MyFile.close();
}



int main(void)
{
    LatticeBoltzmann Ondas;
    int t, tmax=100;
    double rho0=0, Jx0=0, Jy0=0;

    //Start
    Ondas.Star(rho0, Jx0, Jy0);
    //Run
    for  (t = 0; t < tmax; t++)
    {
        Ondas.Collision();
        Ondas.ImposeFields(t);
        Ondas.Adveccion();
    }
    //show
    Ondas.Print("ondas.dat");

    return 0;
}


