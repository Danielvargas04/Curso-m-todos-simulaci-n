#include <iostream>
#include <fstream>
#include <cmath>

//Constantes del problemas

#define Lx 128
#define Ly 128
#define Q 5
#define N 32 //Threads per Block
const int M=(Lx*Ly+N-1)/N; //Blocks per Grid

const int ArraySize=Lx*Ly*Q;

const double W0=1.0/3;

const double C=0.5; //C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
//------ PROGRAMMING ON THE DIVICE (GPU)--------

__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];

__constant__ float d_C[3]; //C[i]: C, C2, AUX0
__constant__ float d_tau[3]; //tau[i]: tau, Utau, UmUtau

//--------Funtions by device--------------

__device__ int d_n(int ix, int iy, int i)
{
    return (ix*Ly+iy)*Q+i;
}

__device__ float d_rho(int ix, int iy, float *d_f)
{
    float sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=d_n(ix,iy,i);
        sum+=d_f[n0];
    }
    return sum; 
}


__device__ float d_Jx(int ix, int iy, float *d_f)
{
    float sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=d_n(ix,iy,i);
        sum+=d_Vx[i]*d_f[n0];
    }
    return sum; 
}

__device__ float d_Jy(int ix, int iy, float *d_f)
{
    float sum=0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=d_n(ix,iy,i);
        sum+=d_Vy[i]*d_f[n0];
    }
    return sum; 
}

__device__ float d_feq(float rho0, float Jx0, float Jy0, int i)
{
    return 3*d_w[i]*(d_C[1]*rho0+d_Vx[i]*Jx0+d_Vy[i]*Jy0);
}

__device__ float d_f0eq(float rho0, float Jx0, float Jy0)
{
    return rho0*d_C[2];
}

//-------------------------Kernels------------------------
__global__ void d_Collision(float *d_f, float *d_fnew)
{
    int ix, iy, i, n0, icell;
    float rho0, Jx0, Jy0;
    //Find which thread and which cell should I work
    icell=blockIdx.x*blockDim.x+threadIdx.x;
    ix=icell/Ly; iy=icell%Ly;
    //Compute the macroscopic fields
    rho0 = d_rho(ix, iy, d_f);
    Jx0 = d_Jx(ix, iy, d_f);
    Jy0 = d_Jy(ix, iy, d_f);
    //Collide and compute fnew
    n0=d_n(ix,iy,0);
    d_fnew[n0]=d_tau[2]*d_f[n0] + d_tau[1]*d_f0eq(rho0, Jx0, Jy0);
    for (i = 1; i < Q; i++)
    {
        n0=d_n(ix,iy,i);
        d_fnew[n0]=d_tau[2]*d_f[n0] + d_tau[1]*d_feq(rho0, Jx0, Jy0, i);
    }
}

__global__ void d_ImposeFields(float *d_f, float *d_fnew, float RhoSource)
{
    int ix, iy, i, n0;
    float rho0, Jx0, Jy0;
    //Find which thread and which cell should I work
    ix=Lx/2; iy=Ly/2;
    //Compute the macroscopic fields
    rho0 = RhoSource;
    Jx0 = d_Jx(ix, iy, d_f);
    Jy0 = d_Jy(ix, iy, d_f);
    //Collide and compute fnew
    n0=d_n(ix,iy,0);
    d_fnew[n0]=d_f0eq(rho0, Jx0, Jy0);
    for (i = 1; i < Q; i++)
    {
        n0=d_n(ix,iy,i);
        d_fnew[n0]=d_feq(rho0, Jx0, Jy0, i);
    }
}

__global__ void d_Advection(float *d_f, float *d_fnew)
{
    int ix, iy, i, n0, icell, ixnew, iynew, n0new;
    //Find which thread and which cell should I work
    icell = blockIdx.x*blockDim.x+threadIdx.x;
    ix=icell/Ly; iy=icell%Ly;
    //Move the contents to the neighboring cells
    for (i = 0; i < Q; i++)
    {
        ixnew=(ix+d_Vx[i]+Lx)%Lx; iynew=(iy+d_Vy[i]+Ly)%Ly;
        n0=d_n(ix,iy,i); n0new=d_n(ixnew,iynew,i);
        d_f[n0new]=d_fnew[n0];
    }
}


//---------------class Laticce--------------
class LatticeBoltzmann
{
private:
    float h_C[3]; //C[i]: C, C2, AUX0
    float h_tau[3]; //tau[i]: tau, Utau, UmUtau
    float h_w[Q];           //Weights
    int h_Vx[Q], h_Vy[Q];    //Velocity vectors
    float *h_f, *h_fnew;    
    float *d_f, *d_fnew;    
public:
    LatticeBoltzmann(void);
    ~LatticeBoltzmann();
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    //------------ Campos macroscopicos-----------------
    float h_rho(int ix, int iy);
    //------------ Funciones de equilibrio-----------------
    float feq(float rho0, float Jx0, float Jy0, int i);
    //---------- Evolucion temporal----------------
    void Star(float rho0, float Jx0, float Jy0);
    void Collision(void);
    void ImposeFields(int t);
    void Adveccion(void);
    //---------- Funciones Globales----------------
    void Print(const  char * NameFile);
};

//-------------Implementacion de funciones-----------------

LatticeBoltzmann::LatticeBoltzmann(void)
{
    //Set constans
    h_C[0]=C; h_C[1]=C2; h_C[2]=AUX0;
    h_tau[0]=tau; h_tau[1]=Utau; h_tau[2]=UmUtau; 
    //Set the weights
    h_w[0]=W0; h_w[1]=h_w[2]=h_w[3]=h_w[4]=(1.0-W0)/4;
    //Set Velocity vectors
    h_Vx[0]=0; h_Vx[1]=1; h_Vx[2]=0; h_Vx[3]=-1; h_Vx[4]=0;
    h_Vy[0]=0; h_Vy[1]=0; h_Vy[2]=1; h_Vy[3]=0; h_Vy[4]=-1;
    //------ Send to the divice------
    cudaMemcpyToSymbol(d_w, h_w,    Q*sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_Vx, h_Vx,  Q*sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_Vy, h_Vy,  Q*sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_C, h_C,    3*sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_tau, h_tau,3*sizeof(float), 0, cudaMemcpyHostToDevice);

    //Create the dynamic arrays
    h_f=new float [ArraySize]; h_fnew=new float [ArraySize];
    //Build the dynamic matrices on the device
    cudaMalloc((void**) &d_f, ArraySize*sizeof(float));
    cudaMalloc((void**) &d_fnew, ArraySize*sizeof(float));
    
}

LatticeBoltzmann::~LatticeBoltzmann()
{
    delete[] h_f; delete[] h_fnew;
    cudaFree(d_f); cudaFree(d_fnew);
}

float LatticeBoltzmann::h_rho(int ix, int iy)
{
    float sum=0.0;
    int n0;
    for (int i = 0; i < Q; i++)
    {
        n0=n(ix,iy,i);
        sum+=h_fnew[n0];

    }
    return sum;    
}

float LatticeBoltzmann::feq(float rho0, float Jx0, float Jy0, int i)
{
    if(i>0)
        return 3*h_w[i]*(C2*rho0 + h_Vx[i]*Jx0 + h_Vy[i]*Jy0);
    else
        return rho0*AUX0;
}

void LatticeBoltzmann::Star(float rho0, float Jx0, float Jy0)
{
    int ix, iy, i, n0;
    for(ix = 0; ix < Lx; ix++)      //for each cell
        for(iy = 0; iy < Ly; iy++)
            for(i = 0; i < Q; i++)  //on each direction
            {
                n0=n(ix,iy,i);
                h_f[n0]=feq(rho0, Jx0, Jy0, i);
            }
    //send to the device
    cudaMemcpy(d_f, h_f, ArraySize*sizeof(float), cudaMemcpyHostToDevice);
}

void LatticeBoltzmann::Collision(void)
{
    //Do by device

    dim3 ThreadsPerBlock(N,1,1);
    dim3 BlocksPerGrid(M,1,1);
    d_Collision<<<BlocksPerGrid, ThreadsPerBlock>>>(d_f, d_fnew);
}

void LatticeBoltzmann::ImposeFields(int t)
{
    //Do by device  
    float lambda=10, omega=2*M_PI/lambda*C;
    float RhoSource=10*sin(omega*t);
    dim3 ThreadsPerBlock(1,1,1);
    dim3 BlocksPerGrid(1,1,1);
    d_ImposeFields<<<BlocksPerGrid, ThreadsPerBlock>>>(d_f, d_fnew, RhoSource);
}

void LatticeBoltzmann::Adveccion(void)
{
    //Do by device

    dim3 ThreadsPerBlock(N,1,1);
    dim3 BlocksPerGrid(M,1,1);
    d_Advection<<<BlocksPerGrid, ThreadsPerBlock>>>(d_f, d_fnew);
}

//--------------Funciones Globales--------------

void LatticeBoltzmann::Print(const char * NameFile)
{
    std::ofstream MyFile(NameFile); double rho0; int ix,iy;
    //Bring back the data from device to host
    cudaMemcpy(h_fnew, d_fnew, ArraySize*sizeof(float), cudaMemcpyDeviceToHost);
    for (ix = 0; ix < Lx; ix++)
    {
        for (iy = 0; iy < Ly; iy++)
        {
            rho0=h_rho(ix, iy);
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


