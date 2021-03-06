#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

typedef struct df_transport{

	int nx;
	double Lx;
	
	double dx, dt, tmax, tnow;
	
	double *xx;
	
	double *un;
	double *unp1;
	double *uex;
	
} df_transport;

extern df_transport null_dft;

void init_df_transport(df_transport *dft, int nx, double Lx);

void exact_sol(double x, double t, double *u);

void plot_data(df_transport *dft);

void compute_df(df_transport *dft, double tmax, double beta);

void compute_exact(df_transport *dft, double tmax);


double f(double x);

double df(double x);

double fnum(double a, double b);

#endif

