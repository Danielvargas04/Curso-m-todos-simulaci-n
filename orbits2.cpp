#include <iostream>
#include <cmath>
#include<fstream>

//constantes
const double g=9.8; 
const double GM=1.0; 

//Declaracion de la clase
class objet;

class objet
{
private:
    double x, y, vx, vy, fx, fy, mass, radius; 
public:
    void init(double x0, double y0, double vx0, double vy0,
     double mass0, double radius0);
    void calcularfuerza(void);
    void move(double dt);
    double getx(void){return x;};
    double gety(void){return y;};
};


//Implementacion de las funciones
void objet::init(double x0, double y0, double vx0, double vy0, 
    double mass0, double radius0)
{
    x=x0; y=y0; vx=vx0; vy=vy0; mass=mass0; radius=radius0;
}
void objet::calcularfuerza(void)
{
    double aux = GM*mass*pow(x*x + y*y,-1.5);
    fx = -aux*x;    fy = -aux*y;

}
void objet::move(double dt)
{
    x = x + vx* dt; y += vy*dt;
    vx +=fx*dt/mass; vy+=fy*dt/mass;
}

int main(){
    std::ofstream archivo("./data/datosparabola2.dat");
    double t, dt=0.01;
    double r0=10.0;
    double omega = sqrt(GM/pow(r0,3));
    double T=2*M_PI/omega;
    double V0=omega*r0;

    objet tierra;  //nombra el objeto antes creado
    tierra.init(r0, 0, 0, V0, 1, 1);

    for ( t = 0; t < 1.2*T; t+=dt)
    {
        archivo<<tierra.getx()<<" "<<tierra.gety()<<std::endl;
        tierra.calcularfuerza();
        tierra.move(dt);
    }
    archivo.close();
    return 0;
}