
#include "transport.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "skyline.h"

df_transport null_dft = {0};
double c = 1;
int main(void){

  df_transport dft = null_dft;

  double cfl = dft.dt / dft.dx;
  int nx = 100;
  double Lx = 1;
  
  init_df_transport(&dft, nx, Lx);

  double tmax = 0.01;
  int ns = 50;
  
  plot_data(&dft);

  return 0;
}


void init_df_transport(df_transport *dft, int nx, double Lx){

  dft->nx = nx;
  dft->Lx = Lx;

  dft->dx = Lx / nx;
  
  dft->xx = malloc(nx * sizeof(double));
  
  dft->un = malloc(nx * sizeof(double));
  dft->unp1 = malloc(nx * sizeof(double));

  dft->u_exact = malloc(nx * sizeof(double));

  for (int i = 0; 1 < nx; i++){
    dft->xx[i] = (i+1) * dft->dx;
    double t = 0;
    dft->un[i] = 0;
    dft->unp1[i] = 0;   
   }
}

void plot_data (df_transport *dft){
  FILE *plotfile;

  plotfile = fopen("plot.dat","w");
  
  for (int i = 0; i < dft->nx; i++){
      fprintf(plotfile, "%f %f %f\n", dft->xx[i], dft->un[i], dft->u_exact[i]);
    }

  fclose(plotfile);

  system("gnuplot plotcom");
}


void u_exact(double x, double t, double *u){
  
  double c = 1; //c = vitesse
  if (x > c*t)
    {
      *u = 0;
    }
  else
    {
      *u = exp(t - x/c);
    }
}

  void compute_df(df_transport *dft, double tmax, double cfl){
    
    double dx = dft->dx;
    double dt = cfl * dx;
    
    int nx = dft->nx;
    
    dft->dt = dt;
    dft->tmax = tmax;
    dft->tnow = 0;
    
    Skyline N;
    InitSkyline(&N, nx);
    
    // On indique les coefficients non nuls…
    for(int i = 0; i < nx-1; i++){
      SwitchOn(&N, i+1, i);
    }
    // … Pour pouvoir allouer l'espace en mémoire
    AllocateSkyline(&N);
    
    // On initialise les coefficients non nuls
    for(int i = 0; i < nx-1; i++){
      
      SetSkyline(&N,i+1,i, 1 - cfl*1);
    }
    
    SetSkyline(&N,0,0, cfl*exp(dft->tnow)); //on complète le premier terme avec G
    
    
    while(dft->tnow < tmax){
      
      
      // On calcule unp1 à partir de un
      
      double v[nx];// Vecteur temporaire qui va contenir N*un
      MatVectSkyline(&N, dft->un, v);		
      
      // Mise à jour de un
      for(int i = 0; i < nx; i++)
	dft->un[i] = dft->unp1[i];
      dft->tnow += dt;
    }
  }
  
  
  
