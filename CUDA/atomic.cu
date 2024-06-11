#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Lx 31
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__global__ void ReduceThreads(float * d_a,float * d_ReduceResult){
  //Set global index
  int ix=threadIdx.x+blockIdx.x*blockDim.x;
  //Define and load shared memory
  extern __shared__ unsigned int temp[];
  temp[threadIdx.x]=d_a[ix];
  __syncthreads();
  //Reduce loop
  for(unsigned int s=blockDim.x/2;s>0;s>>=1){
    if(threadIdx.x<s)
      temp[threadIdx.x]+=temp[threadIdx.x+s];
    __syncthreads();
  }
  //Write to d_ReduceResult
  if(threadIdx.x==0)
  {
    d_ReduceResult[blockIdx.x]=temp[0];
  }
}
int main(){
  //DECLARAR LAS MATRICES
  int ix,bx;
  //DECLARAR LAS MATRICES
  //En el Host
  float h_a[Lx],h_ReduceResult[Mx];
  //En el Device
  float *d_a;             cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float *d_ReduceResult;  cudaMalloc((void**) &d_ReduceResult,Mx*sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++) h_a[ix]=ix;
  for(bx=0;bx<Mx;bx++) h_ReduceResult[bx]=0;
  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_ReduceResult,h_ReduceResult,Mx*sizeof(float),cudaMemcpyHostToDevice);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_a[ix]<<" "; cout<<endl;
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  ReduceThreads<<<BlocksPerGrid,ThreadsPerBlock,Mx*sizeof(float)>>>(d_a,d_ReduceResult);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_ReduceResult,d_ReduceResult,Mx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(bx=0;bx<Mx;bx++) cout<<h_ReduceResult[bx]<<" "; cout<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);  cudaFree(d_ReduceResult);

  return 0;
}