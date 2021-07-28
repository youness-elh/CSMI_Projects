#include "transport.h"
#include "skyline.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

df_transport null_dft = {0};

//Condition au bord gauche
double g(double t)
{
  return exp(-t);
}

int main(void){

	df_transport dft = null_dft;
	
	int nx = 1000;
	double Lx = 1;
	
	init_df_transport(&dft, nx, Lx);
	//printf("%f\n",compute_mass(&dft));
	
	double tmax = 0.6;
	double beta = 0.1;// On prend c=1, donc beta=dt/dx
	compute_df(&dft, tmax, beta);
	compute_exact(&dft, dft.tnow);
	//printf("%f\n",compute_mass(&dft));
	
	plot_data(&dft);

	return 0;
	
}

void init_df_transport(df_transport *dft, int nx, double Lx){

	dft->nx = nx;// (*dft).nx = nx
	dft->Lx = Lx;
	
	dft->dx = Lx / nx;
	
	dft->xx = malloc(nx * sizeof(double));
	
	dft->un = malloc(nx * sizeof(double));
	dft->unp1 = malloc(nx * sizeof(double));
	dft->uex = malloc(nx * sizeof(double));
	
	for(int i = 0; i < nx; i++){
		dft->xx[i] = (i + 1) * dft->dx;
		exact_sol(dft->xx[i], 0, &(dft->un[i]));
		dft->unp1[i] = dft->un[i];
		dft->uex[i] = 0;
	}
}

void exact_sol(double x, double t, double *u){

	if (x >= t)
		*u = 0;
	else
		*u = g(t-x);
}


void plot_data(df_transport *dft){

	double err_sup = 0;
	double err_L2 = 0;
	double err_L1 = 0;
	double errloc = 0;

	FILE *plotfile;
	
	plotfile = fopen("plot.dat", "w");
	
	fprintf(plotfile, "x un us");
	for(int i = 0; i < dft->nx; i++){
		fprintf(plotfile, "%f %f %f\n", dft->xx[i], dft->un[i], dft->uex[i]);
		errloc= fabs(dft->uex[i] - dft->un[i]);
		err_sup = errloc > err_sup ? errloc : err_sup;
		err_L2 += dft->dx * errloc * errloc;
		err_L1 += dft->dx * errloc;
	}
	err_L2 = sqrt(err_L2);
	printf("erreur L2=%f L1=%f\n", err_L2, err_L1);
	
	fclose(plotfile);
	
	system("gnuplot plotcom");
}


void compute_df(df_transport *dft, double tmax, double beta){

	double dx = dft->dx;
	double dt = beta * dx;
	
	int nx = dft->nx;
	
	dft->dt = dt;
	dft->tmax = tmax;
	dft->tnow = 0;
	
	while(dft->tnow < tmax){
	  dft->unp1[0] = exp(- dft->tnow);
	  for(int i = 1; i < nx; ++i){
	    dft->unp1[i] = (1 - beta)*dft->un[i] + beta*dft->un[i-1];
	  }
	  for(int i=0; i<nx; ++i){
	    dft->un[i] = dft->unp1[i];
	  }
	  dft->tnow += dt;
	}
}

void compute_exact(df_transport *dft, double tmax)
{
	for (int i = 0; i < dft->nx; i++){
		exact_sol(dft->xx[i], tmax, &(dft->uex[i]));
	}
}
