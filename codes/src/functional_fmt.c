#include "../include/functional_fmt.h"

#include<stdio.h>
#include<math.h>
#include <string.h>
#include <stdlib.h>

#include "../include/globals.h"
#include "../include/init.h"
#include "../include/utils.h"


fftw_complex  *omega3=NULL;
fftw_complex  *omega2=NULL;
fftw_complex  *omega1=NULL;
fftw_complex  *omega0=NULL;

fftw_complex *omega1v=NULL;
fftw_complex *omega2v=NULL;

double *n0=NULL;
double *n1=NULL;
double *n2=NULL;
double *n3=NULL;
double *n1v=NULL;
double *n2v=NULL;

double *dphidn0=NULL;
double *dphidn1=NULL;
double *dphidn2=NULL;
double *dphidn3=NULL;
double *dphidn1v=NULL;
double *dphidn2v=NULL;

void get_eta()
{
  eta=0.;
  for(int i=0;i<Nspecies;i++)
    eta+=(PI/6.)*rhob[i];
}

void get_mu_HS()
{
  double contribution_HS=(14.*eta-13.*eta*eta+5.*eta*eta*eta)/(2.*(1.-eta)*(1.-eta)*(1.-eta)) - log(1.-eta);
  for(int i=0;i<Nspecies;i++)
    mu[i]+= log(rhob[i]) + contribution_HS;
  
}

void initialize_FMT_omegas()
{
  omega3=fftw_alloc_complex((size_t)K);
  omega2=fftw_alloc_complex((size_t)K);
  omega1=fftw_alloc_complex((size_t)K);
  omega0=fftw_alloc_complex((size_t)K);
  
  omega1v=fftw_alloc_complex((size_t)K);
  omega2v=fftw_alloc_complex((size_t)K);
  
  
  //initialize all omega(k)
  double dk=2.*PI/Lx;
	double k=0.;
	omega3[0][0]=4./3*PI*R*R*R;
	omega3[0][1]=0.;
	
	omega2[0][0]=4.*PI*R*R;
	omega2[0][1]=0.;
	
	omega2v[0][0]=0.0;
	omega2v[0][1]=0.;
	
	for(int i=1;i<K;i++)
	{
		k=dk*i;
		omega2[i][0]=4.*PI*R/k*sin(k*R);
		omega2[i][1]=0.;
		
		omega3[i][0]=4.*PI/(k*k*k)*(sin(k*R)-k*R*cos(k*R));
		omega3[i][1]=0.;
		
		omega2v[i][0]=0.0;
		omega2v[i][1]=(-4.*PI)*(-k*R*cos(k*R)+sin(k*R))/(k*k);
	}
	
	
	for(int i=0;i<K;i++)
	{
		omega1[i][0]=omega2[i][0]/(4.*PI*R);
		omega1[i][1]=0.0;
		
		omega0[i][0]=omega2[i][0]/(4.*PI*R*R);
		omega0[i][1]=0.0;
		
		omega1v[i][0]=0.0;
		omega1v[i][1]=omega2v[i][1]/(4.*PI*R);
	}
}

void initialize_FMT_n()
{
  n0 = malloc(N * sizeof(double));
  n1 = malloc(N * sizeof(double));
  n2 = malloc(N * sizeof(double));
  n3 = malloc(N * sizeof(double));
  
  n1v = malloc(N * sizeof(double));
  n2v = malloc(N * sizeof(double));
  
  for(int i=0;i<N;i++)
	{
		n0[i]=0.;
		n1[i]=0.;
		n2[i]=0.;
		n3[i]=0.;
		n1v[i]=0.;
		n2v[i]=0.;
	}
}

void initialize_FMT_dphidn()
{
  dphidn0 = malloc(N * sizeof(double));
  dphidn1 = malloc(N * sizeof(double));
  dphidn2 = malloc(N * sizeof(double));
  dphidn3 = malloc(N * sizeof(double));
  dphidn1v = malloc(N * sizeof(double));
  dphidn2v = malloc(N * sizeof(double));
  
  for(int i=0;i<N;i++)
	{
		dphidn0[i]=0.;
		dphidn1[i]=0.;
		dphidn2[i]=0.;
		dphidn3[i]=0.;
		dphidn1v[i]=0.;
		dphidn2v[i]=0.;
	}
}

void getn()
{
	conv_FFT(rhocopy, omega0, n0);
	conv_FFT(rhocopy, omega1, n1);
	conv_FFT(rhocopy, omega2,n2);
	conv_FFT(rhocopy, omega3, n3);
	conv_FFT(rhocopy, omega1v, n1v);
	conv_FFT(rhocopy, omega2v, n2v);
}

void getc1_fmt()
{
	
		for(int i=0;i<N;i++)
			{
				dphidn0[i] = -log(1-n3[i]);
				dphidn1[i] =  n2[i]/(1-n3[i]);
				dphidn2[i] = n1[i]/(1-n3[i]) + (3*n2[i]*n2[i] - 3*(n2v[i]*n2v[i]))/(24*PI*(1-n3[i])*(1-n3[i]));
				dphidn3[i] = n0[i]/(1-n3[i]) + (n1[i]*n2[i] - n1v[i]*n2v[i])/(1-n3[i])/(1-n3[i]) + (n2[i]*n2[i]*n2[i] - 3*n2[i]*(n2v[i]*n2v[i]))/12/PI/(1-n3[i])/(1-n3[i])/(1-n3[i]);
				dphidn1v[i] = -n2v[i]/(1-n3[i]);
				dphidn2v[i] = -n1v[i]/(1-n3[i]) - 3*n2[i]*n2v[i]/(12*PI*(1-n3[i])*(1-n3[i]));
			}
	
	
	conv_FFT(dphidn0, omega0, c1_temp);	
	for(int i=0;i<Nspecies;i++)
	  add_c(&(c1[IDX(i,0)]),c1_temp);
	
	conv_FFT(dphidn1, omega1, c1_temp);	
	for(int i=0;i<Nspecies;i++)
	  add_c(&(c1[IDX(i,0)]),c1_temp);
	  
	conv_FFT(dphidn2, omega2, c1_temp);	
	for(int i=0;i<Nspecies;i++)
	  add_c(&(c1[IDX(i,0)]),c1_temp);
	
  conv_FFT(dphidn3, omega3, c1_temp); 
	for(int i=0;i<Nspecies;i++)
	  add_c(&(c1[IDX(i,0)]),c1_temp);
	
	conv_FFT(dphidn1v, omega1v, c1_temp);	
	for(int i=0;i<Nspecies;i++)
	  add_c(&(c1[IDX(i,0)]),c1_temp);
	
	conv_FFT(dphidn2v, omega2v, c1_temp);	
	for(int i=0;i<Nspecies;i++)
	  add_c(&(c1[IDX(i,0)]),c1_temp);
	

}

