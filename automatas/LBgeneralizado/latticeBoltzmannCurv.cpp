#ifndef __latticeBoltzmannCurv__CPP__INCLUDED
#define __latticeBoltzmannCurv__CPP__INCLUDED
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

double Automata::Fuerza(int coord, int ix, int iy, int iz, int t){
      double suma=0.0;
      double delta;
        if(ix!=0&&ix!=Lx-1&&iy!=0&&iy!=Ly-1&&iz!=0&&iz!=Lz-1)
        suma=4.0*Derivada(coord,0,ix,iy,iz,t)+4.0*Derivada(coord,1,ix,iy,iz,t)+4.0*Derivada(coord,2,ix,iy,iz,t);
   
       return suma/2.0;
  /*double suma=0.0;
   double delta;
   if(ix==0||ix==Lx-1){
     if(iy==0||iy==Ly-1){
       if(iz==0||iz==Lz-1)
	 {
	   //es Vertice
	 }
       else
	 {
	   //es Arista z
	   suma=4.0*Derivada(coord,2,ix,iy,iz,t);
	 }
     }
     else{
       if(iz==0||iz==Lz-1)
	 {
	   //Es Arista y
	   suma=4.0*Derivada(coord,1,ix,iy,iz,t);
	 }
       else
	 {
	   //es Cara x
	   suma=4.0*Derivada(coord,1,ix,iy,iz,t)+4.0*Derivada(coord,2,ix,iy,iz,t);
	 }
     }
   }
   else{
     if(iy==0||iy==Ly-1){
       if(iz==0||iz==Lz-1)
	 {
	   //arista x
	   suma=4.0*Derivada(coord,0,ix,iy,iz,t);
	 }
       else
	 {
	   //cara y
	   suma=4.0*Derivada(coord,0,ix,iy,iz,t)+4.0*Derivada(coord,2,ix,iy,iz,t);
	 }
     }
     else{
       if(iz==0||iz==Lz-1)
	 {
	   //Cara z
	   suma=4.0*Derivada(coord,0,ix,iy,iz,t)+4.0*Derivada(coord,1,ix,iy,iz,t);
	 }
       else
	 {
	   //Punto interior
	   suma=4.0*Derivada(coord,0,ix,iy,iz,t)+4.0*Derivada(coord,1,ix,iy,iz,t)+4.0*Derivada(coord,2,ix,iy,iz,t);
	 }
     }
   }
   return suma/2.0;*/
}

double Automata::Derivada(int coord, int j, int ir,int it,int iz, int t){
  
  double delta;
  if(coord==j) delta=1; else delta=0;
  if(j==0){
    //                                                    i      j   r          W
    return (sdetg[ir+1][it][iz]*rho(ir+1, it, iz, t)*(cs2*delta-c2*ginv[coord][j][ir+1][it][iz])/8.0
	    -sdetg[ir-1][it][iz]*rho(ir-1, it, iz,t)*(cs2*delta-c2*ginv[coord][j][ir-1][it][iz])/8.0); 
    //Los dos términos se deben a la suma sobre n, solo hay dos velocidades en cada direccion.
  }
  else if(j==1)
    {
      //                                                          i      j   r        W
      return (sdetg[ir][it+1][iz]*rho(ir, it+1, iz,t)*(cs2*delta-c2*ginv[coord][j][ir][it+1][iz])/8.0
	      -sdetg[ir][it-1][iz]*rho(ir, it-1, iz,t)*(cs2*delta-c2*ginv[coord][j][ir][it-1][iz])/8.0); 
      //Los dos términos se deben a la suma sobre n, solo hay dos velocidades en cada direccion.
    }
  else
    {
      //                                                    i      j   r          W
      return (sdetg[ir][it][iz+1]*rho(ir, it, iz+1,t)*(cs2*delta-c2*ginv[coord][j][ir][it][iz+1])/8.0
	      -sdetg[ir][it][iz-1]*rho(ir, it, iz-1,t)*(cs2*delta-c2*ginv[coord][j][ir][it][iz-1])/8.0); 
      //Los dos términos se deben a la suma sobre n, solo hay dos velocidades en cada direccion.
    }
}



double Automata::Uz(int ix,int iy, int iz, int t){
  int i; double suma;
  for(suma=0,i=0;i<q;i++)
    suma+=V[2][i]*f[actual][ix][iy][iz][i];
  return (suma+Fuerza(2, ix,iy, iz,t)-(c2*rho(ix,iy,iz,t)/2.0)*(FacP[ix][iy][iz]))/sdetg[ix][iy][iz];
}
double Automata::Uy(int ix,int iy, int iz, int t){
  int i; double suma;
  for(suma=0,i=0;i<q;i++)
    suma+=V[1][i]*f[actual][ix][iy][iz][i];
  return (suma+Fuerza(1, ix,iy, iz,t)-(c2*rho(ix,iy,iz,t)/2.0)*(FacT[ix][iy][iz]))/sdetg[ix][iy][iz];
}
double Automata::Ux(int ix,int iy, int iz, int t){
  int i; double suma;
  for(suma=0,i=0;i<q;i++)
    suma+=V[0][i]*f[actual][ix][iy][iz][i];
  return (suma+Fuerza(0, ix,iy, iz,t)-(c2*rho(ix,iy,iz,t)/2.0)*(FacR[ix][iy][iz]))/sdetg[ix][iy][iz];
}

double Automata::rho(int ix,int iy, int iz, int t){
  int i; double suma;
  for(suma=0,i=0;i<q;i++)
    suma+=f[actual][ix][iy][iz][i];
  if(iy==1&&(ix<=Lx/2+2&&ix>=Lx/2-2)&&((iz>=2*Lz/3-2&&iz<=2*Lz/3+2))) /*||(iz>=Lz/3-2&&iz<=Lz/3+2))*/
    return (suma/sdetg[ix][iy][iz]+Amp*sin(omega*t));//*(sin(2.0*M_PI*iz/double(Lz-1))));
  else if(iy==1&&(ix<=Lx/2+2&&ix>=Lx/2-2)&&((iz>=Lz/3-2&&iz<=Lz/3+2))&&t>2000)
    return (suma/sdetg[ix][iy][iz]+Amp*sin(omega*(t-2000)+M_PI));//*(sin(2.0*M_PI*iz/double(Lz-1))));
  else
    return suma/sdetg[ix][iy][iz];
  
}

double Automata::feq(double rho0,double Ux0,double Uy0, double Uz0, int i, int ix, int iy, int iz){
 double UdotVi=Ux0*V[0][i]+Uy0*V[1][i]+Uz0*V[2][i];  
   if(i==0)
     return sdetg[ix][iy][iz]*rho0*w[0];
   else
   return w[i]*sdetg[ix][iy][iz]*(UdotVi/cs2+rho0);

}

void Automata::Muestre(int FrameID, int t){
  int ix,iy, iz;

  double r, Thet, Phi;
  
  string Name1="./Data/CocleaEscalada.";
  string Result1;//string which will contain the result
  stringstream convert1; // stringstream used for the conversion
  convert1 << FrameID;//add the value of FrameID to the characters in the stream
  Result1 = convert1.str();//set Result to the content of the stream
  string Ext=".vtk";
  string Frame1=Name1+Result1+Ext;

  string Name2="./Data/MembranBasilar.";
  string Result2;//string which will contain the result
  stringstream convert2; // stringstream used for the conversion
  convert2 << FrameID;//add the value of FrameID to the characters in the stream
  Result2 = convert2.str();//set Result to the content of the stream
  string Frame2=Name2+Result2+Ext;

   string Name3="./DatosParaInterpolacion/Membran.";
  string Result3;//string which will contain the result
  stringstream convert3; // stringstream used for the conversion
  convert3 << FrameID;//add the value of FrameID to the characters in the stream
  Result3 = convert3.str();//set Result to the content of the stream
  string Frame3=Name3+Result3+Ext;

  
  
  ofstream outClientFile, PuntosMembrana, MembranaInterp;
  outClientFile.precision(6);
  PuntosMembrana.precision(6);
  MembranaInterp.precision(6);
 
  outClientFile.open(Frame1.c_str(), ios::out);
  PuntosMembrana.open(Frame2.c_str(), ios::out);
    MembranaInterp.open(Frame3.c_str(), ios::out);

    outClientFile<<"# vtk DataFile Version 2.0"<<endl;
    outClientFile<<"Archivo VTK para graficar las ondas en geometrías generales "<<endl;
    outClientFile<<"ASCII"<<endl;
    outClientFile<<endl;
    outClientFile<<"DATASET STRUCTURED_GRID"<<endl;
    outClientFile<<"DIMENSIONS "<<Lx<<" "<<Ly<<" "<<Lz<<" "<<endl;
    outClientFile<<"POINTS "<<Lx*Ly*Lz<<" float"<<endl;
    outClientFile<<endl;
    PuntosMembrana<<"# vtk DataFile Version 2.0"<<endl;
    PuntosMembrana<<"Archivo VTK para graficar la membrana basilar en geometrías generales "<<endl;
    PuntosMembrana<<"ASCII"<<endl;
    PuntosMembrana<<endl;
    PuntosMembrana<<"DATASET STRUCTURED_GRID"<<endl;
    PuntosMembrana<<"DIMENSIONS "<<Lx<<" "<<Ly<<" "<<1<<" "<<endl;
    PuntosMembrana<<"POINTS "<<Lx*Ly*1<<" float"<<endl;
    PuntosMembrana<<endl;
    for(iz=0;iz<Lz;iz++){
      for(iy=0;iy<Ly;iy++){
	for(ix=0;ix<Lx;ix++){
	  r=rmin+double(ix)*DR;
	  Thet=tmin+double(iy)*DT;
	  Phi=pmin+double(iz)*DP;
	  outClientFile<<fixed<<cos(Thet)*(R-a*Thet+r*cos(Phi)*(1+cos(Phi))*cos(Thet/F2))*rmaxreal/rmax
		       <<"\t "<<sin(Thet)*(R-a*Thet+r*cos(Phi)*(1+cos(Phi))*cos(Thet/F2))*rmaxreal/rmax
		       <<"\t "<<(r*sin(Phi)*(1+cos(Phi))*cos(Thet/F2)+Thet*A)*rmaxreal/rmax<<endl; //enviar estos datos al archivo.
	  if(iz==Lz/2+1){
	    PuntosMembrana<<fixed<<cos(Thet)*(R-a*Thet+r*cos(Phi)*(1+cos(Phi))*cos(Thet/F2))*rmaxreal/rmax
			  <<"\t "<<sin(Thet)*(R-a*Thet+r*cos(Phi)*(1+cos(Phi))*cos(Thet/F2))*rmaxreal/rmax
			  <<"\t "<<(r*sin(Phi)*(1+cos(Phi))*cos(Thet/F2)+Thet*A)*rmaxreal/rmax<<endl; //enviar estos datos al archivo.
	    MembranaInterp<<fixed<<r*rmaxreal/rmax
			  <<"\t "<<Thet
			  <<"\t "<<rho(ix,iy,iz,t)+rho(ix,iy,iz-2,t)<<endl; //enviar estos datos al archivo.
	  }
	}
      }
      
    }
    outClientFile<<"POINT_DATA "<<Lx*Ly*Lz<<endl;
    outClientFile<<endl;
    outClientFile<<"SCALARS Pressure float 1"<<endl;
    outClientFile<<"LOOKUP_TABLE default"<<endl;
    PuntosMembrana<<"POINT_DATA "<<Lx*Ly*1<<endl;
    PuntosMembrana<<endl;
    PuntosMembrana<<"SCALARS Pressure float 1"<<endl;
    PuntosMembrana<<"LOOKUP_TABLE default"<<endl;
    for(iz=0;iz<Lz;iz++){
      for(iy=0;iy<Ly;iy++){
	for(ix=0;ix<Lx;ix++){
	  r=rmin+double(ix)*DR;
	  Thet=tmin+double(iy)*DT;
	  outClientFile<<fixed<<rho(ix,iy,iz,t)<<endl;  //enviar estos datos al archivo.
	  if(iz==Lz/2)
	    //if(r<=(rmin+double(Lx/2)*DR)+(13.0*(rmax-rmin)/(36.0*5.0*M_PI)*Thet+(rmax-rmin)/18.0)&&r>=(rmin+double(Lx/2)*DR)-(13.0*(rmax-rmin)/(36.0*5.0*M_PI)*Thet+(rmax-rmin)/18.0))
	      PuntosMembrana<<fixed<<rho(ix,iy,iz+1,t)+rho(ix,iy,iz-1,t)<<endl;
	  //else
	  // PuntosMembrana<<fixed<<10000<<endl;
	}
      }
    }
    
    outClientFile.close();
    PuntosMembrana.close();
    MembranaInterp.close();
}
void Automata::Inicie(){
  int i,ix,iy,iz;
  actual=0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	for(i=0;i<q;i++)
	  f[actual][ix][iy][iz][i]=feq(0,0,0,0,i,ix,iy,iz);
}
void Automata::Evolucione(double t){
  double rho0,Ux0,Uy0,Uz0;
  int ix,iy,iz,i, siguiente;
  int ixn, iyn, izn;
  int ixb, iyb, izb;
  siguiente=(actual+1)%2;
  double Sigma=3.0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){

	//Calculo las variables macroscopicas
	rho0=rho(ix,iy,iz,t);  Ux0=Ux(ix,iy,iz,t);  Uy0=Uy(ix,iy,iz,t); Uz0=Uz(ix,iy,iz,t);
	
	for(i=0;i<q;i++){
	
	  ixn=ix+V[0][i];
	  iyn=iy+V[1][i];
	  izn=iz+V[2][i];
	  
	  if(ixn>=0&&ixn<=Lx-1&&izn>=0&&izn<=Lz-1&&iyn>=0&&iyn<=Ly-1&&(izn>=Lz/2+1||izn<=Lz/2-1))
	    {
	      f[siguiente][ixn][iyn][izn][i]=f[actual][ix][iy][iz][i]
		-1.0/Tau*(f[actual][ix][iy][iz][i]-feq(rho0,Ux0,Uy0,Uz0,i,ix, iy, iz));
	    }
	  // else 
	  else {
	    //  if(ixn>=0&&ixn<=Lx-1&&izn>=0&&izn<=Lz-1&&iyn>=Ly-Ly/20&&iyn<=Ly-1){
	    //   f[siguiente][ixn][iyn][izn][i]=f[actual][ix][iy][iz][i]
	    //     -1.0/Tau*(f[actual][ix][iy][iz][i]-feq(rho0,Ux0,Uy0,Uz0,i,ix, iy, iz));
	    // }
	    if(ix==0)
	      f[siguiente][0][iy][iz][i]=0.5*f[actual][1][iy][iz][i];
	    if(iy==0)
	      f[siguiente][ix][0][iz][i]=0.5*f[actual][ix][1][iz][i];
	    if(iz==0)
	      f[siguiente][ix][iy][0][i]=0.5*f[actual][ix][iy][1][i];
	    
	     if(iz==Lz/2-1)
	       f[siguiente][ix][iy][Lz/2-1][i]=0.5*f[actual][ix][iy][Lz/2-2][i];
	     if(iz==Lz/2+1)
	       f[siguiente][ix][iy][Lz/2+1][i]=0.5*f[actual][ix][iy][Lz/2+2][i];
	    
	    if(ix==Lx-1)
	      f[siguiente][Lx-1][iy][iz][i]=0.5*f[actual][Lx-2][iy][iz][i];
	    if(iy==Ly-1)
	      f[siguiente][ix][Ly-1][iz][i]=0.5*f[actual][ix][Ly-2][iz][i];
	    if(iz==Lz-1)
	      f[siguiente][ix][iy][Lz-1][i]=0.5*f[actual][ix][iy][Lz-2][i];	 	
	    // if(iy==Ly-1&&iz<Lz/2)
	    //   f[siguiente][ix][iy][iz][i]=f[actual][ix][iy][Lz-1-iz][i];
	    
	  }
	}
	
      }
  
  
  actual=siguiente;
}




Automata::Automata(void){
  int ix,iy,iz,i,j,k;
  double r,Phi,Thet;

   //Las Funciones
 
  f=new double**** [2];
  for(j=0;j<2;j++) {
    f[j]=new double*** [Lx];
    for(ix=0;ix<Lx;ix++){
      f[j][ix]=new double** [Ly];
      for(iy=0;iy<Ly;iy++){
	f[j][ix][iy]=new double* [Lz];
	for(iz=0;iz<Lz;iz++){
	  f[j][ix][iy][iz]=new double [q];
	}
      }
    }   
  }

  ginv=new double**** [3];
    for(i=0;i<3;i++){
      ginv[i]=new double*** [3];
      for(j=0;j<3;j++){
	ginv[i][j]=new double** [Lx];
	for(ix=0;ix<Lx;ix++){
	  ginv[i][j][ix]=new double* [Ly];
	  for(iy=0;iy<Ly;iy++){
	    ginv[i][j][ix][iy]=new double [Lz];
	  }
	}
      }
    }
    
    Chr=new double** [3];
	for(i=0;i<3;i++){
	  Chr[i]=new double* [3];
	  for(j=0;j<3;j++){
	    Chr[i][j]=new double [3];
	  }
	}
	
	FacR=new double** [Lx]; FacT=new double** [Lx]; FacP=new double** [Lx];
	for(ix=0;ix<Lx;ix++){
	  FacR[ix]=new double* [Ly]; FacT[ix]=new double* [Ly]; FacP[ix]=new double* [Ly];
	  for(iy=0;iy<Ly;iy++){
	    FacR[ix][iy]=new double [Lz]; FacT[ix][iy]=new double [Lz]; FacP[ix][iy]=new double [Lz];
	  }
	}
	
	
  //Creamos las constantes geometricas
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
      r=rmin+double(ix)*DR;
      Thet=tmin+double(iy)*DT;
      Phi=pmin+double(iz)*DP;
      //sdetg[ix][iy]=(cos(Thet)-1)*(cos(Thet)-1)*R*DR*DT*DZ;
      //Trigonometricas
      double c1=cos(Phi),  c22=cos(Thet/F2), s1=sin(Phi),  s22=sin(Thet/F2);
      //Auxiliares:
      
      double alpha=(R-a*Thet+r*c1*(1.0+c1)*c22);
      double alpha2=alpha*alpha;
      double F22=F2*F2;
      double delta=c1+0.5;
      double bet=1.0+c1;
      double line=(-a*Thet+R);

      //Raiz cuadrada de g
      sdetg[ix][iy][iz]=alpha*r*bet*bet*c22*c22*DP*DR*DT;
      
      //Tensor metrico inverso
      ginv[0][0][ix][iy][iz]=(1/(DR*DR))    *(2.0*r*r*(bet*bet)*(F2*F2*c1*c1-0.5*c1-0.5)*c22*c22+4.0*r*F2*F2*c1*bet*line*c22-4.0*(-a*c1*c1+(A*s1-0.5*a)*c1+0.5*A*s1+0.5*a)*F2*r*bet*s22+
					      (r*r+4.0*(a*a-A*A)*F2*F2)*c1*c1*c1+(-8.0*A*F2*F2*s1*a+3.0*r*r)*c1*c1
					      +3.0*(r*r+(A*A-a*a)*F2*F2)*c1+2.0*A*F2*F2*s1*a+r*r+((2.0*Thet*Thet+1.0)*a*a-4.0*R*a*Thet+A*A+2.0*R*R)*F2*F2)/(c22*c22*F2*F2*bet*bet*bet*alpha2);
      ginv[1][1][ix][iy][iz]=(1/(DT*DT))    *1.0/alpha2;
      ginv[2][2][ix][iy][iz]=(1/(DP*DP))    *(r*r*c1*c1*bet*bet*c22*c22+2.0*r*c1*bet*line*c22+(A*A-a*a)*c1*c1+2.0*A*c1*s1*a+(Thet*Thet+1.0)*a*a-2.0*R*a*Thet+R*R)/(c22*c22*bet*bet*r*r*alpha2);
      
      
      
      ginv[0][1][ix][iy][iz]=(1/(DR*DT))    *(r*bet*bet*s22-2.0*(-a*c1*c1+(A*s1-0.5*a)*c1+0.5*A*s1+0.5*a)*F2)/(c22*F2*bet*bet*alpha2);
      ginv[1][0][ix][iy][iz]=ginv[0][1][ix][iy][iz];
      
      
      ginv[0][2][ix][iy][iz]=(1/(DR*DP))    *(r*r*F2*c1*c1*s1*bet*bet*c22*c22+2*r*F2*c1*s1*bet*line*c22-r*bet*bet*(A*c1+s1*a)*s22+2.0
					      *(-2.0*A*a*c1*c1*c1+((A*A-a*a)*s1-A*a)*c1*c1+((0.5*A*A-0.5*a*a)*s1+1.5*A*a)*c1
						+((0.5*Thet*Thet+0.5)*a*a-R*a*Thet+0.5*R*R)*s1+A*a)*F2)/(c22*c22*F2*r*bet*bet*bet*alpha2);
      ginv[2][0][ix][iy][iz]=ginv[0][2][ix][iy][iz];
      
      
      ginv[1][2][ix][iy][iz]=(1/(DT*DP))    *(-A*c1-a*s1)/(c22*bet*r*alpha2);
      ginv[2][1][ix][iy][iz]=ginv[1][2][ix][iy][iz];
      
      //Simbolos de Christoffel
      Chr[0][0][1]=(DR*DT/DR)   *(-2.0*c1*F2*(-a*c1*c1+(A*s1-0.5*a)*c1+0.5*A*s1+0.5*a)*c22-s22*bet*line)/(c22*alpha*bet*F2);
      Chr[0][1][0]=              Chr[0][0][1];
      Chr[0][1][1]=(DT*DT/DR)   *(-2.0*(F22*c1*c1-0.5*F22*c1-0.5)*bet*bet*bet*r*r*c1*c22*c22-4.0*line*(F22*c1*c1-F22*c1*0.5+0.25)*r*bet*bet*c22
				  +4.0*(-a*c1*c1*c1+(A*s1-a*0.5)*c1*c1+A*c1*s1*0.5-0.5*a)*F2*r*bet*s22-2.0*c1*c1*c1*c1*r*r-6.0*c1*c1*c1*r*r+(-6.0*r*r-2.0*F22*((Thet*Thet+2)*a*a-2*R*a*Thet+R*R))*c1*c1
				  +(4.0*A*F22*s1*a-2.0*r*r-F22*((Thet*Thet+2.0)*a*a-2.0*R*a*Thet+R*R))*c1+2.0*F22*(A*s1*a+(Thet*Thet*0.5+1.0)*a*a-R*a*Thet+R*R*0.5))/(c22*bet*bet*alpha*F22);
      Chr[0][1][2]=(DP*DT/DR)   *(-4.0*r*delta*(0.5*r*s1*bet*s22+(A*c1*c1+(s1*a-0.5*A)*c1-0.5*s1*a-0.5*A)*F2))/(bet*alpha*F2);
      Chr[0][2][1]=              Chr[0][1][2];
      Chr[0][2][2]=(DP*DP/DR)   *(-3.0*r/bet);
      
      Chr[1][0][1]=(DT*DR/DT)   *(bet*c1*c22/alpha);
      Chr[1][1][0]=              Chr[1][0][1];
      Chr[1][1][1]=(DT*DT/DT)   *(-2.0*r*c1*bet*s22-2.0*a*F2)/(alpha*F2);
      Chr[1][1][2]=(DT*DP/DT)   *(-r*(2.0*c1+1.0)*s1*c22/alpha);
      Chr[1][2][1]=              Chr[1][1][2];
      
      
      Chr[2][0][1]=(DR*DT/DP)   *(-c1*(A*c1+a*s1))/(r*alpha);
      Chr[2][1][0]=             Chr[2][1][0];
      Chr[2][0][2]=(DR*DP/DP)   *(1/r);
      Chr[2][2][0]=(DR*DP/DP)   *(1/r);
      Chr[2][1][1]=(DT*DT/DP)   *(r*r*F22*c1*c1*s1*bet*bet*c22*c22+2.0*r*F2*c1*s1*bet*line*c22+2.0*2*c1*bet*(A*c1+a*s1)*s22+2.0*(c1*A*a+0.5*s1*((Thet*Thet+2.0)*a*a-2.0*R*a*Thet+R*R))
				  *F2)/(c22*bet*r*alpha*F2);
      Chr[2][1][2]=(DT*DP/DP)   *((-r*c1*bet*bet*s22+2.0*delta*(A*c1*s1-a*c1*c1+a)*F2)*c22-s22*bet*line)/(c22*bet*alpha*F2);
      Chr[2][2][1]=              Chr[2][1][2];
      Chr[2][2][2]=(DP*DP/DP)   *(-2.0*s1/bet);

      FacR[ix][iy][iz]=2.0*Chr[0][0][1]*ginv[0][1][ix][iy][iz]+2.0*Chr[0][2][1]*ginv[2][1][ix][iy][iz]+Chr[0][1][1]*ginv[1][1][ix][iy][iz]+Chr[0][2][2]*ginv[2][2][ix][iy][iz];

       FacT[ix][iy][iz]=2.0*Chr[1][0][1]*ginv[0][1][ix][iy][iz]+2.0*Chr[1][2][1]*ginv[2][1][ix][iy][iz]+Chr[1][1][1]*ginv[1][1][ix][iy][iz];

       FacP[ix][iy][iz]=2.0*Chr[2][0][1]*ginv[0][1][ix][iy][iz]+2.0*Chr[2][2][0]*ginv[2][0][ix][iy][iz]+2.0*Chr[2][2][1]*ginv[2][1][ix][iy][iz]+Chr[2][1][1]*ginv[1][1][ix][iy][iz]+Chr[2][2][2]*ginv[2][2][ix][iy][iz];
      

	    

    }

  //  lambda=LT/1.86;  omega=/*0.00114516;*/2*M_PI*c/lambda;
  double lambdaReal=130.8; // en mm
  lambda=(250.0/16.0)*lambdaReal;/*TT*R/10;*/  omega=/*0.00114516;*/2*M_PI*c/lambda;

 //Los pesos
    w[0]=1.0/4.0;
   w[1]=w[2]=w[3]=w[4]=w[5]=w[6]=1.0/8.0;

   //wi[19]= {1.0/3.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

  //Los Vectores V[0][i]=Vix , V[1][i]=Viy
   V[0][0]=0;   V[1][0]=0;   V[2][0]=0;
   V[0][1]=1;   V[1][1]=0;   V[2][1]=0; 
   V[0][2]=0;   V[1][2]=1;   V[2][2]=0;
   V[0][3]=-1;   V[1][3]=0;   V[2][3]=0;
   V[0][4]=0;  V[1][4]=-1;   V[2][4]=0; 
   V[0][5]=0;  V[1][5]=0;   V[2][5]=1; 
   V[0][6]=0;  V[1][6]=0;  V[2][6]=-1; 
  

 
  
}

Automata::~Automata(void){
  int ix,iy,iz,i,j,k;
  for(j=0;j<2;j++){
    for(ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
	for(iz=0;iz<Lz;iz++){
	  delete[] f[j][ix][iy][iz];
	}
	delete[] f[j][ix][iy];
      }
      delete[] f[j][ix]; 
    }
    delete[] f[j]; 
  }
  delete[] f;
 
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      for(ix=0;ix<Lx;ix++){
	for(iy=0;iy<Ly;iy++){
	  delete[] ginv[i][j][ix][iy];
	}
	delete[] ginv[i][j][ix];
      }
      delete[] ginv[i][j];
    }
    delete[] ginv[i];
  }
  delete[] ginv;
  
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      delete[] Chr[i][j];
    }
    delete[] Chr[i];
  }
  delete[] Chr;
}

#endif
