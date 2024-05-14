#include <iostream>
#include <fstream>
#include <cmath>

//Constantes del problemas

const int Lx=256;
const int Ly=64;

const int Q=9;

const double tau=0.55;
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
    double feq(double rho0, double Ux0, double Uy0, int i);
    //---------- Evolucion temporal----------------
    void Start(double rho0, double Ux0, double Uy0);
    void Collision(void);
    void ImposeFields(double Ufan);
    void Adveccion(void);
    //---------- Funciones Globales----------------
    void Print(const  char * NameFile, double Ufan);
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

void LatticeBoltzmann::Collision(void)
{
    int ix, iy, i, n0;
    double rho0, Ux0, Uy0;

    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
        {
            rho0 = rho(ix, iy, false);
            Ux0 = Jx(ix, iy, false)/rho0;
            Uy0 = Jy(ix, iy, false)/rho0;
            for(i = 0; i < Q; i++)  //on each direction
            {
                n0=n(ix,iy,i);
                fnew[n0]=UmUtau*f[n0] + Utau*feq(rho0, Ux0, Uy0, i);
            }
        }
}

void LatticeBoltzmann::ImposeFields(double Ufan)
{
    int ix, iy, i, n0;
    int ixc=Lx/8, iyc=Ly/2, R=Ly/5;
    double rho0, R2=R*R;
    
    //Go through all cell, looking if they are fan or obstacle
    for (ix = 0; ix < Lx; ix++)
        for(iy = 0; iy < Ly; iy++)
        {
            rho0=rho(ix, iy, false);
            //fan
            if (ix==0)
                for ( i = 0; i < Q; i++)
                {
                    n0=n(ix, iy, i);
                    fnew[n0]=feq(rho0,Ufan,0,i);
                }
            //Obstacle
            else if ((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2)
                for ( i = 0; i < Q; i++)
                {
                    n0=n(ix, iy, i);
                    fnew[n0]=feq(rho0,0,0,i);
                }
            //An extra point at one side to break the isotropy
            else if (ix==ixc && iy==iyc+R+1)
                for ( i = 0; i < Q; i++)
                {
                    n0=n(ix, iy, i);
                    fnew[n0]=feq(rho0,0,0,i);
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

void LatticeBoltzmann::Print(const char * NameFile, double Ufan)
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
    std::cout << "plot '-' using 1:2:3:4 with vectors head filled lt 2 title 'Campo vectorial'" << std::endl;
    double rho0, Ux0, Uy0; int ix,iy;
    for (ix = 0; ix < Lx; ix+=4)
    {
        for (iy = 0; iy < Ly; iy+=4)
        {
            rho0=rho(ix, iy, true);
            Ux0=Jx(ix, iy, true)/rho0;
            Uy0=Jy(ix, iy, true)/rho0;
            std::cout<<ix<<" "<<iy<<" "<<Ux0/Ufan*4<<" "<<Uy0/Ufan*4<<std::endl;
        }
        std::cout<<std::endl;
    }
    std::cout << "e" << std::endl; // Indica a Gnuplot que ha terminado de recibir datos para este frame
}



int main(void)
{
    LatticeBoltzmann air;
    int taux=0, t, tmax=10000;
    double rho0=1.0, Ufan0=0.1;

    //Start
    air.Start(rho0, Ufan0, 0);


    std::cout<<"set terminal gif animate "<<std::endl;
    std::cout<<"set output 'campo_vectorial.gif'"<<std::endl;
    std::cout<<"set xrange [0:256]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set yrange [0:64]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set size ratio -1    # Para mantener la proporción"<<std::endl;

    //Run
    for  (t = 0; t < tmax; t++)
    {
        air.Collision();
        air.ImposeFields(Ufan0);
        air.Adveccion();
        if(taux==25)
        {
            air.Print("null.dat", Ufan0);
            taux=0;
        }
        taux++;
    }
    //show
    

    

    return 0;
}