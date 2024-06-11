#include<iostream>
#include<fstream>
#include<cmath>
#include<sstream>
#include<string>
#include<ctime>
#include<stdio.h>
#include"Constantes.h"
#include"latticeBoltzmannCurv.h"
using namespace std;

class Intensidades{
private:
  double **Maximo, **Minimo;
public:
  Intensidades(void);
  ~Intensidades(void);
  void Inicie(void);
  void Actualice(Automata & Ondas, int t);
  double DeMaximo(int ix, int iy);
  double DeIntensidad(int ix, int iy);
};

double Intensidades::DeMaximo(int ix, int iy){
  return Maximo[ix][iy];
}

double Intensidades::DeIntensidad(int ix, int iy){
  return pow(0.5*(Maximo[ix][iy]-Minimo[ix][iy]),2);
}

void Intensidades::Inicie(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      Maximo[ix][iy]=Minimo[ix][iy]=0;
}
void Intensidades::Actualice(Automata & Ondas,int t){
  int ix,iy; double rho0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      rho0=Ondas.rho(ix,iy,3*Lz/4,t);
      if(rho0>Maximo[ix][iy]) Maximo[ix][iy]=rho0;  
      if(rho0<Minimo[ix][iy]) Minimo[ix][iy]=rho0;  
    }
}

Intensidades::Intensidades(void){
  int ix;
  Maximo=new double* [Lx];  Minimo=new double* [Lx];
  for(ix=0;ix<Lx;ix++){
    Maximo[ix]=new double [Ly];     Minimo[ix]=new double [Ly]; 
  }
  
}


Intensidades::~Intensidades(void){
  int ix;
  for(ix=0;ix<Lx;ix++){
    delete[] Maximo[ix];    delete[] Minimo[ix];
  }
  delete[] Maximo;  delete[] Minimo;

}




int main(){
  Automata Coclea;
  int t;
  int tMax=16000, FrameTime=tMax/800;
  int i=0;
  Intensidades Medicion;
  ofstream Grafica("Medidas.dat", ios::out);
  ofstream Grafica2("MedidasPresion.dat", ios::out);
  Coclea.Inicie();
  for(t=0;t<=tMax;t++){
    Coclea.Evolucione(double(t));
    Medicion.Actualice(Coclea, double(t));
    Grafica<<t<<"\t"<<Medicion.DeMaximo(Lx/2,Ly/4)<<"\t"<<Medicion.DeMaximo(Lx/2,2*Ly/4)<<"\t"<<Medicion.DeMaximo(Lx/2,3*Ly/4)<<endl;
    Grafica2<<t<<"\t"<<Coclea.rho(Lx/2,Ly/4,3*Lz/4,t)<<"\t"<<Coclea.rho(Lx/2,2*Ly/4,3*Lz/4,t)<<"\t"<<Coclea.rho(Lx/2,3*Ly/4,3*Lz/4,t)<<endl;
    if(t%FrameTime==0&&t>0*FrameTime){
      Coclea.Muestre(i,t);
      i++;
      cout<<"Time:"<<t<<" Frame: "<<i<<endl;
    }
  }   
  return 0;
}



