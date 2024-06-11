#ifndef __latticeBoltzmannCurv_H_INCLUDED
#define __latticeBoltzmannCurv_H_INCLUDED
class Automata{
private:
  int actual;
  double w[q];
  int V[3][q];
  double *****f;//f[Lx][Ly][9]=f[ix][iy][i]
  //double f[2][Lx][Ly][Lz][19];
  double *****ginv;
  double ***Chr, ***FacR, ***FacT, ***FacP;
 double lambda;  double omega;
  double sdetg[Lx][Ly][Lz];
public:
   Automata(void);
  ~Automata(void);
  void Inicie(void);
  void Muestre(int FrameID, int t);
  void Evolucione(double t);
  double feq(double rho0,double Ux0,double Uy0, double Uz0, int i, int ix, int iy, int iz);
  double rho(int ix,int iy, int iz, int t);
  double Ux(int ix,int iy, int iz, int t);
  double Uy(int ix,int iy, int iz, int t);
  double Uz(int ix,int iy, int iz, int t);
  double Fuerza(int coord, int ix, int iy, int iz, int t);
  double Derivada(int coord,int j, int ir,int it,int iz, int t);
};

#endif
