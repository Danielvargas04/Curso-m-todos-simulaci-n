#include <iostream>
#include <fstream>
#include <cmath>


#define Lx 16
#define Nx 8

const int Mx = (Lx+Nx-1)/Nx;

//--------------------KERNELS----------------
__global__ void AddTwoVectors(float *d_a,float *d_b,float *d_c){
 //Which thread should I do?
  int ix;  ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_c[ix]=d_a[ix]+d_b[ix];}

int main()
{
    //DECLARE
    //Declare arrays in the Host
    float h_a[Lx],h_b[Lx],h_c[Lx];
    //Declare arrays in the Device
    float*d_a; cudaMalloc((void**) &d_a,Lx*sizeof(float));
    float*d_b; cudaMalloc((void**) &d_b,Lx*sizeof(float));
    float*d_c; cudaMalloc((void**) &d_c,Lx*sizeof(float));

    //INPUT DATA
    //Set data in the Host
    for(int ix=0;ix<Lx;ix++)
    {
        h_a[ix]=ix; h_b[ix]=2*ix;
    }

    //Send data to the Device
    cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice);

    //PROCESS
    //Run parallel on the Device
    dim3 ThreadsPerBlock(Nx,1,1);
    dim3 BlocksPerGrid(Mx,1,1);
    AddTwoVectors<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

    //SHOW RESULTS
    //Bring back to the Host
    cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);
    for(int ix=0;ix<Lx;ix++)
        std::cout<<ix<<" "<<h_c[ix]<<std::endl;

    //Free dynamic memory
    cudaFree(d_a);  cudaFree(d_b);  cudaFree(d_c);
}