#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"
#include "Random64.h"

//Constantes del problema físico
double Lx=160, Ly=60; 
const int N=200, Ns=80, Ntot=N+3+Ns;
const double g=9.8, KHertz=1.0e4, Gamma=150, Kcundall=500, mu=0.4;

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
  double theta,omega,tau,I;
public:
  //Metodos de la clase
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0); tau=0;};// Inline
  void SumeFuerza(vector3D dF,double dtau){F+=dF; tau+=dtau;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  friend class Colisionador;
};

class Colisionador{
private:
  //Atributos
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];
public:
  //Metodos  
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo * Grano,double dt, int N_live);
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
			  double & xCundall,double & sold,double dt);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
  theta=theta0; omega=omega0; I=2.0/5.0*m*R*R;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt); theta+=omega*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m); omega+=tau*(coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
    <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano,double dt, int N_live){
  int i,j;
  //Borro las fuerzas de todos los granos
  for(i=0; i<Ntot; i++){
    Grano[i].BorreFuerza();
  }
  //Sumo fuerza de gravedad de particulas vivas
  vector3D Fg;
  for(i=0; i<N_live; i++){
    Fg.load(0,-Grano[i].m*g,0);
    Grano[i].SumeFuerza(Fg,0);
  }
  //Calculo la fuerza de cada particula con el entorno
  for(i=0; i<N_live; i++){
    for(j=i+1; j<N_live; j++)//Recorre los pares de particulas vivas
    {
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
    }

    for (j = N; j < Ntot; j++)//Recorre sobre las fronteras
    {
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
    }
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
				      double & xCundall,double & sold,double dt){
  //Determinar si hay colision
  vector3D r21=Grano2.r-Grano1.r; double d=r21.norm();
  double R1=Grano1.R,R2=Grano2.R;
  double s=(R1+R2)-d;

  if(s>0){ //Si hay colisión
    //Vectores unitarios
    vector3D n=r21*(1.0/d),t,k; t.load(n.y(),-n.x(),0); k.load(0,0,1);

    //Calculo la velocidad de contacto
    vector3D Rw; Rw.load(0,0,R2*Grano2.omega+R1*Grano1.omega);
    vector3D Vc=(Grano2.V-Grano1.V)-(Rw^n);
    double Vcn=Vc*n, Vct=Vc*t;

    //Fn (Hertz-Kuramoto-Kano)
    double m1=Grano1.m, m2=Grano2.m; double m12=m1*m2/(m1+m2);
    double Fn=KHertz*std::pow(s,1.5)-Gamma*std::sqrt(s)*m12*Vcn;
    
    //Calculo la fuerza tangencial (Cundall)
    xCundall+=Vct*dt; double Ft=-Kcundall*xCundall; double Ftmax=mu*std::fabs(Fn);
    if(std::fabs(Ft)>Ftmax) Ft=Ft/std::fabs(Ft)*Ftmax;

    //Calcula y Cargue las fuerzas
    vector3D F1,F2,tau1,tau2;
    F2=n*Fn+t*Ft; tau2=((n*(-R2))^F2); F1=F2*(-1); tau1=((n*R1)^F1);
    Grano2.SumeFuerza(F2,tau2*k);   Grano1.SumeFuerza(F1,tau1*k);    
  }
  if(sold>=0 && s<0) xCundall=0;  
  sold=s;
}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<std::endl; 
  std::cout<<"set output 'radios_rand.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-10:"<<Lx+10<<"]"<<std::endl;
  std::cout<<"set yrange[-10:"<<Ly+10<<"]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;  
}
void InicieCuadro(void){
    std::cout<<"plot 0,0 ";
    std::cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    std::cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    std::cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    std::cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

void InicieImagenFinal(Cuerpo * Grano, double x_left, double y_left, double y_up, double x_up, double x_right, double y_right);

int main(){
  Cuerpo Grano[Ntot];
  Colisionador Hertz;
  Crandom ran64(78);
  int i,ix,iy;
  //Parametros de la simulación
  double m0=1.0; double R0=2.0;
  double theta=0, omega0=0, omegamax=8.0;
  double x0=Lx/2,y0=Ly-2*R0,Vx0=0,Vy0=0;
  int N_live; 
  double t_relax = 2*std::sqrt(Ly/g), t_aux=0;
  //Variables auxiliares para las paredes
  double Rpared=100*Lx, Mpared=100*m0;
  double Rs = Lx/(2*Ns);
  //Variables auxiliares para correr la simulacion
  double Ncuadros=1, t , tdibujo, dt=1e-3,
        tmax=5,Tcuadro=t_relax/(5*Ncuadros);

  //INICIO
  InicieAnimacion();
  //Inicializar las paredes
  //------------------(  x0,    y0,Vx0,Vy0,theta0,omega0,    m0,    R0)
  Grano[N+Ns].Inicie(Lx/2, Ly+Rpared, 0, 0, 0, 0,Mpared,Rpared); //Pared arriba		
  Grano[N+Ns+1].Inicie(Lx+Rpared, Ly/2, 0, 0, 0, 0,Mpared,Rpared); //Pared derecha	
  Grano[N+Ns+2].Inicie( -Rpared, Ly/2, 0, 0, 0, 0,Mpared,Rpared); //Pared izquierda
  
  for(int i=0;i<Ns;i++){
    //Pared abajo
    Grano[N+i].Inicie(Rs*(2*i+1),0,0,0,0,0,Mpared,Rs);
  }

  // Inicializar los granos
  N_live=1;
  double t_live;
  Grano[0].Inicie(x0, y0, Vx0, Vy0, theta, omega0, m0, R0);
  std::ofstream datafile;
  datafile.open("radios_rand.dat");

  for(t=tdibujo=0; t<tmax ;t+=dt, tdibujo+=dt, t_aux+=dt){

    //Inicializa las particulas en el tiempo de relajacion
    if( t_aux > t_relax && N_live<=N-1){
      //Omega (y/o radio) aleatorio  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Linea a modificar para obtener el punto c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      omega0 = omegamax*(2*ran64.r()-1);
      //R0 = 0.8*(2.0 + ran64.r());

      Grano[N_live].Inicie(x0, y0, Vx0, Vy0, theta, omega0, m0, R0);
      ++N_live;
      t_aux = 0;
    }

    //Creacion de los cuadros del GIF
    if(tdibujo>Tcuadro){
      InicieCuadro();
      for (int i = 0; i<N+Ns; i++)Grano[i].Dibujese();    
      TermineCuadro();
      tdibujo=0; 
    }
    //Integracion de movimiento
    for(i=0;i<N_live;i++) Grano[i].Mueva_r(dt,xi);    
    Hertz.CalculeTodasLasFuerzas(Grano,dt,N_live); 
    for(i=0;i<N_live;i++) Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N_live;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt,N_live); 
    for(i=0;i<N_live;i++) Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N_live;i++) Grano[i].Mueva_r(dt,Um2chiplusxi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt,N_live); 
    for(i=0;i<N_live;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N_live;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt,N_live); 
    for(i=0;i<N_live;i++)Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N_live;i++) Grano[i].Mueva_r(dt,xi);
  }

  //Analisis de los angulos del perfil de la pila
  double x_left, y_left, y_up, x_up, x_right, y_right, tang, error_rel;
  x_left=Grano[0].Getx();
  x_right=Grano[0].Getx();
  y_up=Grano[0].Gety();

  for (int i = 1; i < N_live; i++)
  {
    if (x_left>Grano[i].Getx()){x_left=Grano[i].Getx(); y_left=Grano[i].Gety();}
    if (x_right<Grano[i].Getx()){x_right=Grano[i].Getx(); y_right=Grano[i].Gety();}
    if (y_up<Grano[i].Gety()){y_up=Grano[i].Gety(); x_up=Grano[i].Getx();}
  }

  tang = (y_up/std::fabs(Lx/2.0-x_left)+y_up/std::fabs(Lx/2.0-x_right))/2;
  error_rel = std::fabs(tang/mu - 1); 

  InicieImagenFinal(Grano, x_left, y_left, y_up, x_up, x_right, y_right);

  datafile<<"Tan"<<'\t'<<"mu"<<'\t'<<"Error relativo"<<'\t'<<"y_up"<<'\t'<<"x_up"<<'\t'<<"y_left"<<'\t'<<"x_left"<<'\t'<<"y_right"<<'\t'<<"x_right"<<std::endl;
  datafile<<tang<<'\t'<<mu<<'\t'<<error_rel<<'\t'<<y_up<<'\t'<<x_up<<'\t'<<y_left<<'\t'<<x_left<<'\t'<<y_right<<'\t'<<x_right<<std::endl;

  datafile.close();
  return 0;
}

void InicieImagenFinal(Cuerpo * Grano, double x_left, double y_left, double y_up, double x_up, double x_right, double y_right){
  // Cambia el terminal a uno adecuado para imágenes estáticas, por ejemplo PNG
  std::cout << "\n set terminal png" << std::endl;
  std::cout << "set output 'radios_rand.png'" << std::endl;
  // Los demás comandos de configuración pueden permanecer igual
   // Puedes llamar a la función InicieCuadro si es adecuada para el cuadro final
  std::cout << "set arrow from " << x_left << "," << y_left << " to " << x_up << "," << y_up << " nohead lc 'red' dashtype 2" << std::endl;
  std::cout << "set arrow from " << x_up << "," << y_up << " to " << x_right << "," << y_right << " nohead lc 'blue' dashtype 2" << std::endl;
  InicieCuadro();
  for (int i = 0; i<N+Ns; i++)Grano[i].Dibujese();
  TermineCuadro(); // Asegura que se envíen los comandos a gnuplot
}