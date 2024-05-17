#include <iostream>
#include <fstream>
#include <cmath>

//Constantes del problemas

const int Lx=1;
const int Ly=64;

const int Q=9;

const double tau=1.2;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double ThreeUmU2tau=3*(1-1/(2*tau));

//---------------class Laticce--------------
class LatticeBoltzmann
{
private:
    double w[Q];            //Weights
    int Vx[Q], Vy[Q];    //Velocity vectors
    double *f, *fnew;            //Grid
public:
    LatticeBoltzmann(void);
    ~LatticeBoltzmann();
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    //------------ Campos macroscopicos-----------------
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew, double Fx);
    double Jy(int ix, int iy, bool UseNew, double Fy);
    //------------ Funciones de equilibrio-----------------
    double feq(double rho0, double Ux0, double Uy0, int i);
    double Fi(int Ux0, int Uy0, double Fx, double Fy, int i);
    //---------- Evolucion temporal----------------
    void Start(double rho0, double Ux0, double Uy0);
    void Collision(double gx, double gy);
    void ImposeFields(void);
    void Adveccion(void);
    //---------- Funciones Globales----------------
    void Print(const  char * NameFile, double gx, double gy);
};

//-------------Implementacion de funciones-----------------

LatticeBoltzmann::LatticeBoltzmann(void)
{
    //Set the weights
    w[0]=4.0/9; 
    w[1]=w[2]=w[3]=w[4]=1.0/9;
    w[5]=w[6]=w[7]=w[8]=1.0/36;

    //Set Velocity vectors
    Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
    Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;

    Vx[5]=1; Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
    Vy[5]=1; Vy[6]=1; Vy[7]=-1; Vy[8]=-1; 
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

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew, double Fx)
{
    double sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        if (UseNew) sum+=Vx[i]*fnew[n0];
        else sum+=Vx[i]*f[n0];
    }
    return sum+0.5*Fx;    
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew, double Fy)
{
    double sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        if (UseNew) sum+=Vy[i]*fnew[n0];
        else sum+=Vy[i]*f[n0];
    }
    return sum+0.5*Fy;    
}

double LatticeBoltzmann::Fi(int Ux0, int Uy0, double Fx, double Fy, int i)
{
    double UdotVi=Ux0*Vx[i]+Uy0*Vy[i];
    double FdotVi=Fx*Vx[i]+Fy*Vy[i];
    double UdotF=Ux0*Fx+Uy0*Fy;
    return ThreeUmU2tau*w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}

double LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i)
{
    double UdotVi=Ux0*Vx[i] + Uy0*Vy[i] , U2=Ux0*Ux0+Uy0*Uy0;
    return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0)
{
    int ix, iy, i, n0;
    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
            for(i = 0; i < Q; i++)  //on each direction
            {
                n0=n(ix,iy,i);
                f[n0]=feq(rho0, Ux0, Uy0, i);
            }
}

void LatticeBoltzmann::Collision(double gx, double gy)
{
    int ix, iy, i, n0;
    double rho0, Ux0, Uy0, Fx, Fy;

    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
        {
            rho0 = rho(ix, iy, false);
            Fx = gx*rho0; Fy = gy*rho0;
            Ux0 = Jx(ix, iy, false, Fx)/rho0;
            Uy0 = Jy(ix, iy, false, Fy)/rho0;
            for(i = 0; i < Q; i++)  //on each direction
            {
                n0=n(ix,iy,i);
                fnew[n0]=UmUtau*f[n0] + Utau*feq(rho0, Ux0, Uy0, i)+ Fi(Ux0, Uy0, Fx, Fy, i);
            }
        }
}

void LatticeBoltzmann::ImposeFields(void)
{
    int ix, iy, i, n0;
    double rho0;
    
    //lower wall
    iy=0;
    for (ix = 0; ix < Lx; ix++)
    {
        rho0=rho(ix, iy, false);
        for (i = 0; i < Q; i++)
        {
            n0=n(ix,iy,i);
            fnew[n0]=feq(rho0, 0, 0, i);
        }
    }
    //Upper wall
    iy=Ly-1;
    for (ix = 0; ix < Lx; ix++)
    {
        rho0=rho(ix, iy, false);
        for (i = 0; i < Q; i++)
        {
            n0=n(ix,iy,i);
            fnew[n0]=feq(rho0, 0, 0, i);
        }
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

void LatticeBoltzmann::Print(const char * NameFile, double gx, double gy)
{
    /*std::ofstream MyFile(NameFile); 
    double rho0, Ux0, Uy0; int ix,iy;
    for (ix = 0; ix < Lx; ix+=4)
    {
        for (iy = 0; iy < Ly; iy+=4)
        {
            rho0=rho(ix, iy, true);
            Ux0=Jx(ix, iy, true)/rho0;
            Uy0=Jy(ix, iy, true)/rho0;
            MyFile<<ix<<" "<<iy<<" "<<Ux0/Ufan*4<<" "<<Uy0/Ufan*4<<std::endl;
        }
        MyFile<<std::endl;
    }
    MyFile.close();*/
    std::cout << "plot '-' with lines"<<std::endl;
    double rho0, Ux0, Uy0, Fx, Fy; 
    int ix=0,iy;
    
    for (iy = 0; iy < Ly; iy++)
    {
        rho0=rho(ix, iy, true);
        Fx=gx*rho0; Fy=gy*rho0;
        Ux0=Jx(ix, iy, true, Fx)/rho0;
        Uy0=Jy(ix, iy, true, Fy)/rho0;
        std::cout<<iy<<" "<<Ux0<<std::endl;
    }
    //std::cout<<std::endl;
    
    std::cout << "e" << std::endl; // Indica a Gnuplot que ha terminado de recibir datos para este frame
}



int main(void)
{
    LatticeBoltzmann air;
    int taux=0, t, tmax=1000;
    double rho0=1.0, g=0.1;

    //Start
    air.Start(rho0, 0, 0);


    std::cout<<"set terminal gif animate "<<std::endl;
    std::cout<<"set output 'forzamiento.gif'"<<std::endl;
    std::cout<<"unset key"<<std::endl;
    std::cout<<"set xrange [0:64]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set yrange [0:150]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set size ratio -1    # Para mantener la proporción"<<std::endl;

    //Run
    for  (t = 0; t < tmax; t++)
    {
        air.Collision(g, 0);
        air.ImposeFields();
        air.Adveccion();
        if(taux==25)
        {
            air.Print("null.dat", g, 0);
            taux=0;
        }
        taux++;
    }
    //show
    //air.Print("null.dat", g, 0);
    return 0;
}