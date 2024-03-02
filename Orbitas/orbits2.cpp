#include <iostream>
#include <cmath>
#include <vector>

//constantes
const double g=9.8; 
const double GM=1.0; 

//Declaracion de la clase
class planet;

class planet
{
private:
    //  posicion, velocidad, fuerza, masa, tamaño
    std::vector<double> r, v, f;
    double mass, radius; 
public:
    void init(double x0, double y0, double z0, double vx0, double vy0,
        double vz0, double mass0, double radius0);
    //Metodos o Acciones
    void calcularfuerza(void);
    void move(double dt);
    void arranque(double dt);
    // Getters
    double getR(int index) const { return r[index]; }
    double getV(int index) const { return v[index]; }
    double getF(int index) const { return f[index]; }

    double getMass() const { return mass; }
    double getRadius() const { return radius; }

    // Setters
    void setR(int index, double value) { r[index] = value; }
    void setV(int index, double value) { v[index] = value; }
    void setF(int index, double value) { f[index] = value; }
};


//Implementacion de las funciones
void planet::init(double x0, double y0, double z0, double vx0, double vy0, 
    double vz0, double mass0, double radius0)
{
    r.resize(3);
    v.resize(3);
    f.resize(3);
    r = {x0, y0, z0};
    v = {vx0, vy0, vz0};
    mass=mass0; radius=radius0;
}
void planet::calcularfuerza(void)
{   
    //fuerza central gravitacional
    double aux = GM*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-1.5);
    for (int i = 0; i <= 2; i++)
    {
       f[i] = -aux*r[i];
    }
}
void planet::arranque(double dt)
{
    for (int i = 0; i <= 2; i++)
    {
        v[i]-=f[i]*dt/(2*mass);
    }
    
}

void planet::move(double dt)
{
    //EDOs resueltas mediante Euler
    for (int i = 0; i <= 2; i++)
    {
        v[i] += f[i]*dt/mass;
        r[i] += v[i]*dt;
    }
}

int main(){
    //constantes de movimiento
    double t, dt=0.001;
    double r0=10.0;
    double omega = sqrt(GM/(pow(r0,3)));
    double T=2*M_PI/omega;
    double V0=omega*r0;

    planet tierra;  //nombra el objeto antes creado
    tierra.init(r0, 0, 0, 0, V0, 0, 1, 1);
    tierra.calcularfuerza();
    tierra.arranque(dt);
    for ( t = 0; t < 1.1*T; t+=dt)
    {
        std::cout<<tierra.getR(0)<<" "<<tierra.getR(1)<<std::endl;
        tierra.calcularfuerza();
        tierra.move(dt);
    }

    return 0;
}