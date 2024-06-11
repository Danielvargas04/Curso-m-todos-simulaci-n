#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lr=200;
const int Ltheta=200;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);
const double Cs2=1.0/4; //This is the value for D2Q5 and D3Q7

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------- Class Declaration ------------
class LatticeBoltzmann;
class Geometry;

//--------------------- class Geometry ------------
class Geometry{
private:
  double r_min,dr,theta_min,dtheta;
public:
  //i,j=0 is r; i,j=1 is theta
  Geometry(void);
  double r(int ir){return r_min+dr*ir;};
  double theta(int itheta){return theta_min+dtheta*itheta;};
  double gRoot(int ir,int itheta);
  double gTrace(int ir,int itheta);
  double MetricTensor_g(int i,int j,int ir,int itheta);
  double InverseMetricTensor_gm1(int i,int j,int ir,int itheta);
  double CristoffelSymbol_Gamma(int i,int j,int k,int ir,int itheta);
  double x(int ir,int itheta);
  double y(int ir,int itheta);
};
Geometry::Geometry(void){
  r_min=0.1; dr=0.1; theta_min=0; dtheta=2*M_PI/Ltheta;
}
double Geometry::gRoot(int ir,int itheta){// $\sqrt{g}$
  double r0=r(ir), theta0=theta(itheta);
  return r0*dr*dtheta;
}
double Geometry::gTrace(int ir,int itheta){
  int i; double sum;
  for(i=0;i<2;i++)
    sum+=MetricTensor_g(i,i,ir,itheta);
  return sum;
}
double Geometry::MetricTensor_g(int i,int j,int ir,int itheta){
  double r0=r(ir), theta0=theta(itheta), g_ij;
  if(i==0 && j==0)
    g_ij=dr*dr;
  else if (i==1 && j==1)
    g_ij=(r0*r0)*(dtheta*dtheta);
  else
    g_ij=0;
  return g_ij;
}
double Geometry::InverseMetricTensor_gm1(int i,int j,int ir,int itheta){
  double r0=r(ir), theta0=theta(itheta), g_ij;
  if(i==0 && j==0)
    g_ij=1.0/(dr*dr);
  else if (i==1 && j==1)
    g_ij=1.0/((r0*r0)*(dtheta*dtheta));
  else
    g_ij=0;
  return g_ij;
}
double Geometry::CristoffelSymbol_Gamma(int i,int j,int k,int ir,int itheta){
  double r0=r(ir), theta0=theta(itheta), Gamma_ijk; //k=superscipt, i,j=subscript
  if(i==1 && j==1 && k==0)
    Gamma_ijk=-r0*(dtheta*dtheta)/dr;
  else if ((i==0 && j==1 && k==1)||(i==1 && j==0 && k==1))
    Gamma_ijk=dr/r0;
  else
    Gamma_ijk=0;
  return Gamma_ijk;
}
double Geometry::x(int ir,int itheta){
  double r0=r(ir), theta0=theta(itheta);
  return r0*cos(theta0);
}
double Geometry::y(int ir,int itheta){
  double r0=r(ir), theta0=theta(itheta);
  return r0*sin(theta0);
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  double w[Q];
  int Vx[Q],Vy[Q]; 
  double *f, *fnew; // f[ix][itheta][i]
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ir,int itheta,int i){return (ir*Ltheta+itheta)*Q+i;};
  double F(int i,int ir,int itheta,bool UseNew,Geometry & Polar);
  double rhog(int ir,int itheta,bool UseNew);
  double rho(int ir,int itheta,bool UseNew,Geometry & Polar);
  double Jxg(int ir,int itheta,bool UseNew,Geometry & Polar);
  double Jyg(int ir,int itheta,bool UseNew,Geometry & Polar);
  double feq(double rhog0,double Jxg0,double Jyg0,int i,double gTrace0);
  void Start(double rhog0,double Jxg0,double Jyg0,Geometry & Polar);
  void Collision(Geometry & Polar);
  void ImposeFields(int t,Geometry & Polar);
  void Advection(void);
  void Print(const char * NombreArchivo,Geometry & Polar);
};  
LatticeBoltzmann::LatticeBoltzmann(void){
  //set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //creat the dynamic arrays
  int ArraySize=Lr*Ltheta*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
double LatticeBoltzmann::F(int a,int ir,int itheta,bool UseNew,Geometry & Polar){
  double metricforce,divergence,delta_ab,Aab_next; int b,c,i,irnext,ithetanext;
  //THE FORCEMENT HAS TWO TERMS:
  //The first term is metricforce, with metricforce=0.5*C^2\rho \sum_{ij}\Gamma^a_{ij}(g^{-1})^{ij}
  for(metricforce=0,b=0;b<2;b++)
    for(c=0;c<2;c++)
      metricforce+=Polar.CristoffelSymbol_Gamma(b,c,a,ir,itheta)*Polar.InverseMetricTensor_gm1(b,c,ir,itheta);
  metricforce*=C2/2*rho(ir,itheta,UseNew,Polar);
  
  //The second term is 1/2*divergence, with divergence=\nabla\cdot A=\partial_b A^{a,b}
  //, where A^{a,b}=c^2\sqrt{g}\rho [ (g^{-1})^{a,b}-\delta^{a,b} ]
  
  //The divergence on A is computed "a la" lattice-Boltzmann,
  //\partial_b A^{a,b}=\frac{1}{c_s^2}\sum_{b,i} w_i v_{i,b} A^{a,b}(\vec x+\vec v_i)
  for(divergence=0,b=0;b<2;b++){
    if(b==a) delta_ab=1; else delta_ab=0;
    for(i=0;i<Q;i++){
      irnext=(ir+Vx[i]+Lr)%Lr; ithetanext=(itheta+Vy[i]+Ltheta)%Ltheta;
      Aab_next=rhog(irnext,ithetanext,UseNew)*(Cs2*delta_ab-C2*Polar.InverseMetricTensor_gm1(a,b,irnext,ithetanext));
      if(b==0) divergence+=w[i]*Vx[i]*Aab_next; else if(b==1) divergence+=w[i]*Vy[i]*Aab_next;
    }
    divergence*=0.5/Cs2;
  }
  return divergence-metricforce;
}
double LatticeBoltzmann::rhog(int ir,int itheta,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ir,itheta,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
double LatticeBoltzmann::rho(int ir,int itheta,bool UseNew,Geometry & Polar){
  double rhog0=rhog(ir,itheta,UseNew);
  double gRoot0=Polar.gRoot(ir,itheta);
  return rhog0/gRoot0;
}
double LatticeBoltzmann::Jxg(int ir,int itheta,bool UseNew,Geometry & Polar){
  double sum,rho0; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ir,itheta,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum+F(0,ir,itheta,UseNew,Polar);
}  
double LatticeBoltzmann::Jyg(int ir,int itheta,bool UseNew,Geometry & Polar){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ir,itheta,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum+F(1,ir,itheta,UseNew,Polar);
}  
double  LatticeBoltzmann::feq(double rhog0,double Jxg0,double Jyg0,int i,double gTrace0){
 double JgdotVi=Jxg0*Vx[i]+Jyg0*Vy[i];
  if(i==0)
    return rhog0*w[0];
  else
    return w[i]*(rhog0+(JgdotVi/Cs2));
}  
void LatticeBoltzmann::Start(double rhog0,double Jxg0,double Jyg0,Geometry & Polar){
  int ir,itheta,i,n0; double gTrace0;
  for(ir=0;ir<Lr;ir++) //for each cell
    for(itheta=0;itheta<Ltheta;itheta++){
      gTrace0=Polar.gTrace(ir,itheta);
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ir,itheta,i);
	f[n0]=feq(rhog0,Jxg0,Jyg0,i,gTrace0);
      }
    }
}  
void LatticeBoltzmann::Collision(Geometry & Polar){
  int ir,itheta,i,n0; double rhog0,Jxg0,Jyg0,gTrace0;
  for(ir=0;ir<Lr;ir++) //for each cell
    for(itheta=0;itheta<Ltheta;itheta++){
      //compute the macroscopic fields
      rhog0=rhog(ir,itheta,false);Jxg0=Jxg(ir,itheta,false,Polar); Jyg0=Jyg(ir,itheta,false,Polar);
      gTrace0=Polar.gTrace(ir,itheta);
      for(i=0;i<Q;i++){ //en cada dirección
	n0=n(ir,itheta,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(rhog0,Jxg0,Jyg0,i,gTrace0);
      }
    }  
}
void LatticeBoltzmann::ImposeFields(int t,Geometry & Polar){
  int i,ir,itheta,n0,n0next; double rhog0,Jxg0,Jyg0,gTrace0;
  double lambda,omega; lambda=25; omega=2*M_PI/lambda*C;
  //A source at the origin
  ir=0; itheta=Ltheta/2;
  rhog0=10*sin(omega*t); Jxg0=Jxg(ir,itheta,false,Polar); Jyg0=Jyg(ir,itheta,false,Polar);
  gTrace0=Polar.gTrace(ir,itheta);
  for(i=0;i<Q;i++){
    n0=n(ir,itheta,i);
    fnew[n0]=feq(rhog0,Jxg0,Jyg0,i,gTrace0);
  }
  //Free surface at the end
  ir=Lr-2; itheta=Ltheta/2;
  for(i=0;i<Q;i++){
    n0=n(ir,itheta,i);  n0next=n(ir+1,itheta,i);
    fnew[n0next]=f[n0];
  }
}
void LatticeBoltzmann::Advection(void){
  int ir,itheta,i,irnext,ithetanext,n0,n0next;
  for(ir=0;ir<Lr;ir++) //for each cell
    for(itheta=0;itheta<Ltheta;itheta++)
      for(i=0;i<Q;i++){ //on each direction
	irnext=(ir+Vx[i]+Lr)%Lr; ithetanext=(itheta+Vy[i]+Ltheta)%Ltheta;
	n0=n(ir,itheta,i); n0next=n(irnext,ithetanext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
void LatticeBoltzmann::Print(const char * FileName,Geometry & Polar){
  ofstream MyFile(FileName); double rho0; int ir,itheta=Ltheta/2;
  for(ir=0;ir<Lr;ir++){
    rho0=rho(ir,itheta,true,Polar);
    MyFile<<Polar.r(ir)<<" "<<rho0<<endl;
    }
  MyFile.close();
}
//------------------- Global Functions ------------

int main(void){
  LatticeBoltzmann Waves;
  Geometry Polar;
  int t,tmax=200;
  double rhog0=0,Jxg0=0,Jyg0=0;

  //Start
  Waves.Start(rhog0,Jxg0,Jyg0,Polar);
  //Run
  for(t=0;t<tmax;t++){
    Waves.Collision(Polar);
    Waves.ImposeFields(t,Polar);
    Waves.Advection();
  }
  //Show
  Waves.Print("WavesOnPolarCoordinates.dat",Polar);
 
  return 0;
}  
