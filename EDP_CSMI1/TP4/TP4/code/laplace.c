#include "laplace.h"
#include "skyline.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <unistd.h>

//#define FULL

df_laplace null_dfl = {0};

int main(void)
{
  df_laplace dfl = null_dfl;

  int nx = 10; //N
  int ny = 1000; //M

  double Lx = 1;

  init_df_laplace(&dfl, nx, ny, Lx);
  
  struct timeval tv1,tv2;             // Calcul du
  double temps;                       // temps
  gettimeofday(&tv1,NULL);            // d'exécution
   
  compute_df(&dfl);

  gettimeofday(&tv2,NULL);            // de la
  temps=(tv2.tv_sec-tv1.tv_sec);      // résolution numérique
  printf("N = %d, M = %d\n",nx,ny);
  printf("temps=%f secondes\n",temps); 

  
  plot_data(&dfl);
  
  return 0;
}



void init_df_laplace(df_laplace *dfl, int nx, int ny, double Lx)
{
	dfl->nx = nx;
	dfl->ny = ny;
	int N = nx*ny;
	dfl->Lx = Lx;
	
	dfl->h = Lx / (nx + 1);
	
	dfl->u = malloc(N * sizeof(double));
	dfl->F = malloc(N * sizeof(double));
	
	
	for(int k = 0; k < N; ++k){
	  dfl->u[k] = 0;
	  dfl->F[k] = dfl->h * dfl->h;
	  //printf("%f", dfl->F[k]);
	}
}


void plot_data(df_laplace *dfl){

	FILE *plotfile;
	plotfile = fopen("plot.dat", "w");

	int nx = dfl->nx;
	double h = dfl->h;
	
	for( int i = 0; i < dfl->nx; ++i ){
	  for( int j = 0; j < dfl->ny; ++j ){
	    fprintf(plotfile, "%f %f %f\n", (i+1)*h, (j+1)*h, dfl->u[i + j*nx]);
	  }
	}
	
	fclose(plotfile);
	
	system("gnuplot plotcom");
}


void compute_df(df_laplace *dfl){

	double h = dfl->h;
	int nx = dfl->nx;
	int ny = dfl->ny;
	int N = nx*ny;
	Skyline A;
	
	InitSkyline(&A, N);
	
	// On indique les coefficients non nuls…
	for(int i = 0; i < N-1 ; ++i){
	  
	  SwitchOn(&A, i, i+1);
	  SwitchOn(&A, i+1, i);
	  
	  if (i+nx < N){
	    SwitchOn(&A, i, i+nx);
	    SwitchOn(&A, i+nx, i);
	  }
	}
	#ifdef FULL
	for (int i=0; i<N; ++i ){
	  for (int j=0; j<N; ++j ){
	    SwitchOn(&A, i, j);
	  }
	}
	#endif
	// … Pour pouvoir allouer l'espace en mémoire
	AllocateSkyline(&A);
	
	// On initialise les coefficients non nuls
	for(int i = 0; i < N ; ++i){

	  
	  SetSkyline(&A,i,i,4);
	  if ((i-(nx-1))%nx != 0){
	  SetSkyline(&A,i,i+1,-1);	
	  SetSkyline(&A,i+1,i,-1);
	  }
	  if (i+nx < N){
	    SetSkyline(&A,i+nx,i,-1);
	    SetSkyline(&A,i,i+nx,-1);
	  }
	
	}

	// On affiche la matrice A
	
	/*for (int i = 0; i<N; ++i){
	  for (int j = 0; j<N; ++j){
	    printf("%0.0f ", GetSkyline(&A,i,j));
	  }
	  printf("\n");
	  }*/
	//DisplaySkyline(&A);
	
	// On calcule la factorisation LU de la matrice qu'on veut inverser
	FactoLU(&A);
    
	FastSolveSkyline(&A, dfl->F, dfl->u);
	

	
}

