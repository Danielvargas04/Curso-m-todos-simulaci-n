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

const double sigma=Ly/9, mu=Lx/2;
const double D=(1/3.0)*(tau-0.5);
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
    double Var(void);
    //------------ Funciones de equilibrio-----------------
    double feq(double rho0, double Ux0, double Uy0, int i);
    //---------- Evolucion temporal----------------
    void Start(double mu, double sigma, double Ux0, double Uy0);
    void Collision(void);
    void ImposeFields(double Ufan);
    void Adveccion(void);
    //---------- Funciones Globales----------------
    void Print(void);
    void Animate(void);
    void Plot(void);
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
    return 0;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew)
{
    /*double sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        if (UseNew) sum+=Vy[i]*fnew[n0];
        else sum+=Vy[i]*f[n0];
    }
    return sum;*/  
    return 0;
}

double LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i)
{
    double UdotVi=Ux0*Vx[i] + Uy0*Vy[i] , U2=Ux0*Ux0+Uy0*Uy0;
    return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

void LatticeBoltzmann::Start(double mu, double sigma, double Ux0, double Uy0)
{
    double rho0;
    int ix, iy, i, n0;
    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
            for(i = 0; i < Q; i++)  //on each direction
            {
                rho0 = (1/(sigma*pow(2*M_PI, 0.5))) *exp(-0.5*pow((ix-mu)/sigma,2));
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

double LatticeBoltzmann::Var(void)
{
    double N, x_prom, Var;
    double aux_prom, aux_Var;

    for (int ix = 0; ix < Lx; ix++)
    {
        for (int iy = 0; iy < Ly; iy++)
        {
            N+=rho(ix, iy, true);
            aux_prom+=rho(ix, iy, true)*ix;
        }
    }
    x_prom = aux_prom/N;

    for (int ix = 0; ix < Lx; ix++)
    {
        for (int iy = 0; iy < Ly; iy++)
        {
            aux_Var+=rho(ix, iy, true)*pow(ix-x_prom,2);
        }
    }    
    Var = aux_Var/N;
    return Var;
}

//--------------Funciones Globales--------------

void LatticeBoltzmann::Print(void)
{
    double rho0; int ix,iy;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for (ix = 0; ix < Lx; ix+=1)
    {
        for (iy = 0; iy < Ly; iy+=1)
        {
            rho0=rho(ix, iy, true);
            std::cout<<ix<<" "<<iy<<" "<<rho0<<std::endl;
        }
        std::cout<<std::endl;
    }
    std::cout << "e" << std::endl;
}

void LatticeBoltzmann::Animate(void)
{
    std::cout<<"set terminal gif animate "<<std::endl;
    std::cout<<"set output 'densidad.gif'"<<std::endl;
    std::cout<<"set xrange [0:256]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set yrange [0:64]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set size ratio -1    # Para mantener la proporción"<<std::endl;
}

void LatticeBoltzmann::Plot(void)
{
    std::cout<<"set terminal jpeg enhanced"<<std::endl;
    std::cout<<"set output 'grafica_var.jpg'"<<std::endl;
    std::cout<<"set xlabel 'Tiempo'"<<std::endl;
    std::cout<<"set ylabel 'Varianza'"<<std::endl;
    std::cout<<"set title 'Varianza en funcion del tiempo'"<<std::endl;
    std::cout<<"D = "<<D<<std::endl;
    std::cout<<"b = "<<sigma*sigma<<std::endl;
    std::cout<<"plot '-' using 1:2 with points title 'Datos Simulacion', \
     2*D*x + b  with lines title 'Ajuste teórico'"<<std::endl;
    
}

int main(void)
{
    LatticeBoltzmann air;
    int taux=0, t, tmax=1000;
    double varianza;
    

    //Start
    double UX0 =0.0, UY0 =0.0; 
    air.Start(mu, sigma, UX0, UY0);

    air.Plot();

    //Run
    for  (t = 0; t < tmax; t++)
    {
        air.Collision();
        //air.ImposeFields(Ufan0);
        air.Adveccion();

        if(taux==int(1000/10))
        {
            varianza=air.Var();
            std::cout<<t<<" "<<varianza<<std::endl;

            taux=0;
        }
        taux++;
    }
    std::cout<<"e"<<std::endl;

    //show
    //air.Print();
    return 0;
}