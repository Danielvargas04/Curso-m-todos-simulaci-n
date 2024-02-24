#include <iostream>
#include <cmath>

//constantes
const double g=9.8; 
const double GM=1.0; 

//Declaracion de la clase
class planet;

class planet
{
private:
    double r[3], v[3], f[3]; double mass, radius; 
public:
    void init(double x0, double y0, double z0, double vx0, double vy0,
        double vz0, double mass0, double radius0);
    void calcularfuerza(void);
    void move(double dt);
    double getx(void){return r[0];};
    double gety(void){return r[1];};
    double getfx(void){return f[0];};
};


//Implementacion de las funciones
void planet::init(double x0, double y0, double z0, double vx0, double vy0, 
    double vz0, double mass0, double radius0)
{
    r[0]=x0;r[1]=y0;r[2]=z0; 
    v[0]=vx0;v[1]=vy0;v[2]=vz0;
    mass=mass0; radius=radius0;
}
void planet::calcularfuerza(void)
{
    double aux = GM*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-1.5);
    for (int i = 0; i <= 2; i++)
    {
       f[i] = -aux*r[i];
    }
}
void planet::move(double dt)
{
    for (int i = 0; i <= 2; i++)
    {
        v[i] += f[i]*dt/mass;
        r[i] += v[i]*dt;
    }
}

int main(){
    double t, dt=0.001;
    double r0=10.0;
    double omega = sqrt(GM/(pow(r0,3)));
    double T=2*M_PI/omega;
    double V0=omega*r0;

    planet tierra;  //nombra el objeto antes creado
    tierra.init(r0, 0, 0, 0, V0, 0, 1, 1);

    for ( t = 0; t < 1.1*T; t+=dt)
    {
        std::cout<<tierra.getx()<<" "<<tierra.gety()<<std::endl;
        tierra.calcularfuerza();
        tierra.move(dt);
    }

    return 0;
}