#include "../include/functional_ES.h"

#include<stdio.h>
#include<math.h>
#include <string.h>
#include <stdlib.h>


#include "../include/globals.h"
#include "../include/init.h"
#include "../include/utils.h"



double lambdaB=0.;
double Vq_L=0.;
double Vq_R=0.;
double *q=NULL;
double *phi=NULL;
double *psi=NULL;
int N_ES=0;
char* BC=NULL;

void initialize_ES_df()
{
  BC=malloc(5*sizeof(char));
  lambdaB=0.0;
  q=malloc(Nspecies*sizeof(double));
  
  N_ES=(iend-NiR+1);
  
  phi=malloc(N_ES*sizeof(double));
  psi=malloc(N_ES*sizeof(double));
  
  initialize_val(phi,N_ES,0.);
  initialize_val(psi,N_ES,0.);
}

void getc1_ES()
{
  initialize_val(phi,N_ES,0.);
  initialize_val(psi,N_ES,0.);
  
  for(int i=0;i<N_ES;i++)
  {
    for(int j=0;j<Nspecies;j++)
    phi[i]+=-q[j]*rho[IDX(j,i+NiR)];
  }
  
  
  poisson_1D(dx, N_ES,Vq_L,Vq_R,phi,psi,BC);
  
  /*
  FILE *t_w=fopen("../data/psi.dat","w");
  for(int i=0;i<N_ES;i++)
  fprintf(t_w,"%d %f %f\n",i,phi[i],psi[i]);
  fclose(t_w);
  */
  
  
  for(int i=0;i<Nspecies;i++)
  for(int j=NiR;j<=iend;j++)
  {
    c1[IDX(i,j)]+=lambdaB*(-q[i])*psi[j-NiR];
  }
}
