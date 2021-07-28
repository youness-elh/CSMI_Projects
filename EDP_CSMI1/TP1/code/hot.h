#ifndef _HOT_H_
#define _HOT_H_

typedef struct df_chaleur{

	int nx;
	double Lx;
	
	double dx, dt, tmax, tnow;
	
	double *xx;
	
	double *un;
	double *unp1;
	double *us;
	
} df_chaleur;

extern df_chaleur null_dfc;

void init_df_chaleur(df_chaleur *dfc, int nx, double Lx);

void init_data(double x, double *u);

void plot_data(df_chaleur *dfc);

void compute_serie(df_chaleur *dfc, double tmax, int ns);

void compute_df(df_chaleur *dfc, double tmax, double cfl, double theta);

#endif

