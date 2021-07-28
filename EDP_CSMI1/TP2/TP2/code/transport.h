#ifndef _TRANSPORT_H
#define _TRANSPORT_H

typedef struct df_transport{

  int nx;
  double Lx;
  
  double dx, dt, tmax, tnow;

  double *xx;
  
  double *un;
  double *unp1;

  double *u_exact;


} df_transport;

extern df_transport null_dfc;

void init_df_transport(df_transport *dft, int nx, double Lx);

void u_exact(double x, double t, double *u);

void plot_data(df_transport *dft);

void compute_df(df_transport *dft, double tmax, double cfl);

#endif
