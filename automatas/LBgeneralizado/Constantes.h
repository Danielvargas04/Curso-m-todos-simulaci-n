#ifndef __Constantes_H_INCLUDED__
#define __Constantes_H_INCLUDED__
const int Lx=10, Ly=10*Lx,Lz=Lx+10;
const double Tau=0.5;
const double UX0=0.0, UY0=0.0, UZ0=0.0, RHO0=0.0;
const int q=7;

//CARACTERISTICAS GEOMETRICAS DE LA CÓCLEA
const double A=135.0, a=75, F2=18.67, R=1280.0;
const double cs2=1.0/4.0;
const double c=0.5;    const double c2=c*c; const double Amp=5.0;
const double l=0.0;
const double tmin=0.0; //Correccion que debe hacerse en theta para que no divergan las cantidades geométricas
const double pmin=-M_PI+1.35;
const double rmin=4.0*10.240; const double rmax=10.0*rmin;
const double TR=rmax-rmin, TT=5.0*M_PI-tmin-(0.0+tmin)/*Es lo mismo que 2PI-2corrt pero lo dejo así para acordarme*/, TP=M_PI-1.35-(-M_PI+1.35); //Los tamaños en R, Theta y Zeta.
const double LR=Lx,        LT=Ly,       LP=Lz;    //EL número de celdas convertido a Double
const double rmaxreal=0.9;
const double DR=TR/LR, DT=TT/LT, DP=TP/LP; //Los delta en cada direccion.

#endif
