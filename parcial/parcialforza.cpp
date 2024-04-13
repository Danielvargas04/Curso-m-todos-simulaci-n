#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"
#include "Random64.h"
//Constantes globales
const int N=5;
const double k=1.0, l0 = 10.0, Gamma = 0.05;

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
    void Inicie(double x0,double V0,double m0,double R0);
    void BorreFuerza(void){F.load(0,0,0);};// Inline
    void SumeFuerza(vector3D dF){F+=dF;};// Inline
    void Mueva_r(double dt,double coeficiente);
    void Mueva_V(double dt,double coeficiente);
    void Dibujese(void);
    double Getx(void){return r.x();}; // Inline
    double Gety(void){return r.y();}; // Inline
    void forza(double x){r.load(x,0,0);}; // Inline
    friend class Colisionador;
};

class Colisionador{
private:
    //Atributos
public:
    //Metodos  
    void CalculeTodasLasFuerzas(Cuerpo * bola, double dt);
    void CalculeFuerzaEntre(Cuerpo & bola1,Cuerpo & bola2);
};

//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double V0, double m0,double R0){
    r.load(x0,0,0);  V.load(V0,0,0); m=m0; R=R0;
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
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * bolas, double dt){
    int i;
    //Borro las fuerzas de todas las bolas
    for(i=0; i<N; i++){
        bolas[i].BorreFuerza();
    }

    //Calculo fuerza de rozamiento
    vector3D fr;
    for (i=0; i < N-1; i++)
    { 
        double vel = bolas[i].V.norm();
        fr.load(-Gamma*bolas[i].m*vel ,0,0);
        bolas[i].SumeFuerza(fr);
    }
    
    //Calculo la fuerza de cada particula con su vecina siguiente
    for (i = 0; i < N; i++)
    {       
        CalculeFuerzaEntre(bolas[i],bolas[i+1]);
    }
    
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & bola0,Cuerpo & bola1){
    //Determinar las distancias entre las 2 bolas
    vector3D r10=bola1.r-bola0.r;

    //Calcula el desplazamiento
    double d = r10.norm()-l0;
    double fk = 0; 

    //Calcular la fuerza de resorte entre las dos
    vector3D f,f0,f1;

    //fuerza del resorte
    fk=-k*d;
    f.load(fk,0,0);
    //bola0 siente fuerza en sentido positivo, bola1 en negativo
    f0 = (-1)*f; f1 = f;

    bola0.SumeFuerza(f0); bola1.SumeFuerza(f1);
}

//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
    std::cout<<"set terminal gif animate"<<std::endl; 
    std::cout<<"set output 'resortes.gif'"<<std::endl;
    std::cout<<"unset key"<<std::endl;
    std::cout<<"set xrange[-5:45]"<<std::endl;
    std::cout<<"set yrange[-8:8]"<<std::endl;
    std::cout<<"set size ratio -1"<<std::endl;
    std::cout<<"set parametric"<<std::endl;
    std::cout<<"set trange [0:7]"<<std::endl;
    std::cout<<"set isosamples 12"<<std::endl;  
}
void InicieCuadro(void){
    std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    std::cout<<std::endl;
}


int main()
{
    //Cosntantes de la simulacion
    double m = 1.0, r = 1.0;
    double x00=0, v0=0;

    double x0=0, x4, x1=1*0, x2=std::sqrt(2)*0, x3=1*0;

    double omega = 0.5, A = 0.1;
    //Variables auxiliares para correr la simulacion
    int i;
    double Ncuadros=100, t , tdibujo, dt=1e-2,
        tmax=30,Tcuadro=tmax/Ncuadros;
    //Objetos
    Cuerpo bolas[N];
    Colisionador resorte;

    //Inicializacion

    bolas[0].Inicie(x00, v0, m, r);
    bolas[1].Inicie(x00 + 1*l0 + x1, v0, m, r);
    bolas[2].Inicie(x00 + 2*l0 + x2, v0, m, r);
    bolas[3].Inicie(x00 + 3*l0 + x3, v0, m, r);    
    bolas[4].Inicie(x00 + 4*l0, v0, m, r);

    //Simulacion    
    InicieAnimacion();

    std::ofstream archivo;
    archivo.open("forzado.dat");
    archivo<<"t\tx"<<std::endl;

    for (t=tdibujo=0; t<tmax ; t+=dt, tdibujo+=dt)
    {     
        archivo<<t<<"\t"<<bolas[1].Getx()<<std::endl;


        //Creacion de los cuadros del GIF
        if(tdibujo>Tcuadro){
            InicieCuadro();
            for (int i = 0; i<N; i++)bolas[i].Dibujese();    
            TermineCuadro();
            tdibujo=0; 
        }

        //Calculo del movimiento forzado
        x4 = 40 + A*sin(omega*t);
        bolas[4].forza(x4);

        //Integracion de movimiento para las 3 particulas moviles
        for(i=1; i<N-1; i++) bolas[i].Mueva_r(dt,xi);    
        resorte.CalculeTodasLasFuerzas(bolas, dt); 
        for(i=1; i<N-1; i++) bolas[i].Mueva_V(dt,Um2lambdau2);
        for(i=1; i<N-1; i++) bolas[i].Mueva_r(dt,chi);
        resorte.CalculeTodasLasFuerzas(bolas, dt); 
        for(i=1; i<N-1; i++) bolas[i].Mueva_V(dt,lambda);
        for(i=1; i<N-1; i++) bolas[i].Mueva_r(dt,Um2chiplusxi);
        resorte.CalculeTodasLasFuerzas(bolas, dt); 
        for(i=1; i<N-1; i++)bolas[i].Mueva_V(dt,lambda);
        for(i=1; i<N-1; i++) bolas[i].Mueva_r(dt,chi);
        resorte.CalculeTodasLasFuerzas(bolas, dt); 
        for(i=1; i<N-1; i++)bolas[i].Mueva_V(dt,Um2lambdau2);
        for(i=1; i<N-1; i++) bolas[i].Mueva_r(dt,xi);
    }
    archivo.close();
    return 0;
}
