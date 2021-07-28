#ifndef __LAPLACE_H__
#define __LAPLACE_H__


typedef struct df_laplace
{
  int nx, ny;
  
  double Lx;
  double h;
  
  double *u;
  double *F;
  
} df_laplace;


void init_df_laplace(df_laplace *dfl, int nx, int ny, double Lx);

void compute_df(df_laplace *dfl);

void plot_data(df_laplace *dfl);

#endif
