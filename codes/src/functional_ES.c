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

double *E_exec=NULL;// = d/dx psi(x)
double *P=NULL; //Polarization density
double *Pnew=NULL;
double p=0.; //dipole moment

void initialize_ES_df()
{
  BC=malloc(5*sizeof(char));
  lambdaB=0.0;
  q=malloc(Nspecies*sizeof(double));
  
  N_ES=iend+1;//(iend-NiR+1);
  
  phi=malloc(N_ES*sizeof(double));
  psi=malloc(N_ES*sizeof(double));
  
  if(Is_polar)
  {
    E_exec=malloc(N*sizeof(double));
    initialize_val(E_exec,N,0.); 
  }
  
  initialize_val(phi,N_ES,0.);
  initialize_val(psi,N_ES,0.);
}

void getc1_ES()
{
  initialize_val(phi,N_ES,0.);
  initialize_val(psi,N_ES,0.);
  if(Is_polar)
  initialize_val(E_exec,N,0.); 
  
  for(int i=0;i<N_ES;i++)
  {
    for(int j=0;j<Nspecies;j++)
    phi[i]+=q[j]*lambdaB*rho[IDX(j,i)];//phi(x) : free charge density (=\sum_{i} q_i rho_i(x))
  }
  
  if(strcmp(BC, "NN") != 0)
  poisson_1D(dx, N_ES,Vq_L,Vq_R,phi,psi,BC);
  else
  {
    if(Is_polar)
    psi_calculator_P(dx,N_ES,Vq_L,Vq_R,phi,psi,E_exec,P,lambdaB);
    else
    psi_calculator(dx,N_ES,Vq_L,Vq_R,phi,psi);
  }
  
  /*
  FILE *t_w=fopen("../data/psi.dat","w");
  for(int i=0;i<N_ES;i++)
  fprintf(t_w,"%d %f %f\n",i,phi[i],psi[i]);
  fclose(t_w);
  */
  
  
  for(int i=0;i<Nspecies;i++)
  for(int j=NiR;j<=iend;j++)
  {
    c1[IDX(i,j)]+=(q[i])*psi[j];
  }
}
