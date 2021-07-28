#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "transport.h"

//#define BURG

df_transport null_dft = {0};

int main(void)
{   
  df_transport dft = null_dft;
  
  int nx = 15;
    
  #ifdef BURG
  double Lx = 3;
  #else
  double Lx = 1;
  #endif
  
  init_df_transport(&dft, nx, Lx);
  double tmax = 0.5;
  double beta = 0.5;
  
  compute_df(&dft, tmax, beta);
  #ifdef BURG
  compute_exact(&dft, dft.tnow);
  #endif
  
  plot_data(&dft);
  
  return 0;
  
}

void init_df_transport(df_transport *dft, int nx, double Lx)
{
  dft->nx = nx;// (*dft).nx = nx
  dft->Lx = Lx;
  dft->dx = Lx / nx ;
  
  dft->xx = malloc(nx * sizeof(double));
  dft->un = malloc(nx * sizeof(double));
  dft->unp1 = malloc(nx * sizeof(double));
  dft->uex = malloc(nx * sizeof(double));
  
  for(int i = 0; i < nx; ++i){
    
    #ifdef BURG
    dft->xx[i] = (i + 1) * dft->dx - 1; 
    exact_sol(dft->xx[i], 0, &(dft->un[i]));
    dft->uex[i] = 0;
    #else
    dft->xx[i] = (i + 1) * dft->dx;
    dft->un[i] = 200; //pmax
    dft->un[nx-1] = 0;
    #endif
    
    dft->unp1[i] = dft->un[i];
  }
}

void exact_sol(double x, double t, double *w) //Solution exacte trouvée avec eq.Carac. [Burger]
{
  if (t < 1){
    if (x < t){
      *w = 1;
    }
    else if (x > 1){
      *w = 0;
    }
    else {
      *w = (1-x)/(1-t);
    }
  }
  else{
    if (x < t/2 + 0.5){
      *w = 1;
    }
    else{
      *w = 0;
    }
  } 
}


void plot_data(df_transport *dft)
{
	double err_sup = 0;
	double err_L2 = 0;
	double errloc = 0;

	FILE* plotfile;
	
	#ifdef BURG
	plotfile = fopen("plot_Burg.dat", "w");
	fprintf(plotfile, "x un uex\n");
	#else
	plotfile = fopen("plot_Traf.dat", "w");
       	fprintf(plotfile, "x un\n");
        #endif
	
	for(int i = 0; i < dft->nx; ++i){
	  #ifdef BURG
	  fprintf(plotfile, "%f %f %f\n", dft->xx[i], dft->un[i], dft->uex[i]);
	  errloc= fabs(dft->uex[i] - dft->un[i]);
	  err_sup = errloc > err_sup ? errloc : err_sup;
	  err_L2 += dft->dx * errloc * errloc;
	  #else
	  fprintf(plotfile, "%f %f\n", dft->xx[i], dft->un[i]);
	  #endif
	}
	
	#ifdef BURG
	err_L2 = sqrt(err_L2);
	printf("erreurs : L2=%f L_infini=%f\n", err_L2, err_sup);
	#endif
	
	fclose(plotfile);
	
	#ifdef BURG
	system("gnuplot plotcom_Burg");
	#else
	system("gnuplot plotcom_Traf");
	#endif
}


void compute_df(df_transport *dft, double tmax, double beta)
{
	double dx = dft->dx;
	int nx = dft->nx;
	
	#ifdef BURG
	double dt = beta * dx; //max(f'[w_0(x)]) = 1
	#else
	double dt = beta * dx / 200; //pmax = 200 km/h
	#endif
	
	double left = 0;
	double right = 0;
	
	while (dft->tnow < tmax){
	  #ifdef BURG
	  exact_sol(-1,dft->tnow,&left);
	  exact_sol(nx,dft->tnow,&right);
	  #else
	  left = 200;
	  right = 0;
	  #endif
	  
	  dft->unp1[0] -= (dt/dx)*(fnum(dft->un[0],dft->un[1])-fnum(left,dft->un[0]));
	  
	  for (int i=1; i<nx; ++i){
	    dft->unp1[i] -= (dt/dx)*(fnum(dft->un[i],dft->un[i+1])-fnum(dft->un[i-1],dft->un[i]));
	  }

	  dft->unp1[nx] -= (dt/dx)*(fnum(dft->un[nx],right)-fnum(dft->un[nx-1],dft->un[nx]));
	    
	  for (int i=0; i<nx+1; ++i){
	    dft->un[i] = dft->unp1[i];
	  }
	  dft->tnow += dt;
	}
	
}

void compute_exact(df_transport *dft, double tmax)
{
	for (int i = 0; i < dft->nx; ++i){
		exact_sol(dft->xx[i], tmax, &(dft->uex[i]));
	}
}


#ifdef BURG
double f(double w) //burger
{
  return pow(w,2)/2;
}

double df(double w)
{
  return w;
}

#else
double f(double p) //traffic
{
  double vmax = 130;
  double pmax = 200;
  
  return p*vmax*(1 - p/pmax);
}

double df(double p)
{
  double vmax = 130;
  double pmax = 200;

  return (1 - 2*p/pmax)*vmax;
}
#endif

double fnum(double a, double b)
{
  return 0.5*( f(a)+f(b) ) -0.5*fmax(fabs(df(a)), fabs(df(b)))*(b - a);
}
