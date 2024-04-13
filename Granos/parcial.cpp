#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"
#include "Random64.h"

//Constantes del problema físico
double Lx=10; 
const int N=1, Ns=1, Ntot=N+Ns;
const double g=9.8, KHertz=1.0e4, Gamma=0;

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

//--------------- Declarar las clases-----------
class Cuerpo;
class Colisionador;

//--------- Declarar las interfases de las clases---------
class Cuerpo{
private:
  //Atributos de la clase
  vector3D r,V,F; 
  double m,R; 
public:
  //Metodos de la clase
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  void Sety(double y){r.load(0,y,0);}; // Inline
  friend class Colisionador;
};

class Colisionador{
private:

public:
  //Metodos  
  void CalculeTodasLasFuerzas(Cuerpo * Grano,double dt);
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2, double dt);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------- Funciones de la clase Colisionador --------

void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano,double dt){
  int i,j;
  //Borro la fuerza sobre la pelota
  Grano[0].BorreFuerza();

  //Sumo fuerza de gravedad 
  vector3D Fg;
  Fg.load(0,-Grano[0].m*g,0);
  Grano[0].SumeFuerza(Fg);

  //Calculo la fuerza de la particula con la raqueta
  CalculeFuerzaEntre(Grano[0],Grano[1], dt);

}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2, double dt){
  //Determinar si hay colision
  vector3D r21=Grano2.r-Grano1.r; double d=r21.norm();
  double R1=Grano1.R,R2=Grano2.R;
  double s=(R1+R2)-d;

  if(s>0){ //Si hay colisión
  std::cout<<"#### s: "<<s<<std::endl;
    //Calculo la velocidad de contacto
    vector3D Vc=(Grano2.V-Grano1.V);
    double Vcn=std::fabs(Vc.y());

    //Fn (Hertz-Kuramoto-Kano)
    double m1=Grano1.m;
    double Fn=KHertz*std::pow(s,1.5)-Gamma*std::sqrt(s)*m1*Vcn;
    if (Fn<0){Fn=0;}
    
    //Calcula y Cargue las fuerzas
    vector3D F1;
    F1.load(0,Fn,0);
    Grano1.SumeFuerza(F1);   
  }

}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<std::endl; 
  std::cout<<"set output 'rebote.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-10:10]"<<std::endl;
  std::cout<<"set yrange[-10:60]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;  
}
void InicieCuadro(double y=0){
    std::cout<<"plot "<<-Lx/7.0<<"*t,"<<y<<"";
    std::cout<<" , "<<Lx/7.0<<"*t,"<<y<<"";        //pared de abajo
}
void TermineCuadro(void){
    std::cout<<std::endl;
}


int main(){
  Cuerpo Grano[Ntot];
  Colisionador Hertz;
  Crandom ran64(78);
  int i;
  //Parametros de la simulación
  double m0=1.0; double R0=2.0;
  double x0=0,y0=30,Vx0=0,Vy0=0;


  //Variables auxiliares para las paredes
  double Rpared=10000.0, Mpared=1000*m0;
  double ypared = 0;

  //Variables auxiliares para correr la simulacion
  double Ncuadros=60, t , tdibujo, dt=1e-3,
        tmax=20 ,Tcuadro=tmax/Ncuadros;

  //INICIO
  InicieAnimacion();
  std::ofstream datafile;
  datafile.open("altura.dat");
  datafile <<"t\ty"<< std::endl;

  //Inicializar la pared
  //------------------(  x0,    y0,Vx0,Vy0,    m0,    R0)
  Grano[1].Inicie(0, ypared-Rpared, 0, 0, Mpared, Rpared); //Pared Abajo		
  
  // Inicializar los granos

  Grano[0].Inicie(x0, y0, Vx0, Vy0, m0, R0);

  for(t=tdibujo=0; t<tmax ;t+=dt, tdibujo+=dt){   
    //Creacion de la figura 

    datafile <<t<<"\t"<<Grano[0].Gety()<< std::endl;

    //Creacion de los cuadros del GIF
    if(tdibujo>Tcuadro){
      InicieCuadro();
      Grano[0].Dibujese();    
      TermineCuadro();
      tdibujo=0; 
    }
    //Integracion de movimiento
    Grano[0].Mueva_r(dt,xi);    
    Hertz.CalculeTodasLasFuerzas(Grano,dt); 
    Grano[0].Mueva_V(dt,Um2lambdau2);
    Grano[0].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); 
    Grano[0].Mueva_V(dt,lambda);
    Grano[0].Mueva_r(dt,Um2chiplusxi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); 
    Grano[0].Mueva_V(dt,lambda);
    Grano[0].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); 
    Grano[0].Mueva_V(dt,Um2lambdau2);
    Grano[0].Mueva_r(dt,xi);
  }
  datafile.close();
  return 0;
}

void data_imagen(std::ofstream data, double t, double y){
  data <<t<<","<<y<< std::endl;
}
