#include "../include/functional_LJ.h"

#include<stdio.h>
#include<math.h>
#include <string.h>
#include <stdlib.h>


#include "../include/globals.h"
#include "../include/init.h"
#include "../include/utils.h"


double *eps=NULL;
double *Ufilter=NULL;
fftw_complex *UfilterFFT=NULL;
double rc=0.0;
double rmin=0.0;
int NiLJ=0.0;
double alphaLJ=0.;

void initialize_LJ_df()
{
  eps=malloc(Nspecies*Nspecies*sizeof(double));
  Ufilter=malloc(N*sizeof(double));
  UfilterFFT=fftw_alloc_complex((size_t)K);
  
  
  
  rc=2.5*sigma;
  rmin=pow(2.,1./6)*sigma;
  NiLJ=(int)(rc/dx);
  iend=N-3*NiLJ;
  alphaLJ=(1.171861897)*(4.*PI);
  
  initialize_val(Ufilter,N,0.);
  initialize_val(eps,Nspecies*Nspecies,0.);
  getUfilterFFT();
  
}

double getU(double z)
{
	if(z>=rc)
	return 0.;
	if(z>=rmin)
	return (.4*pow(z,-10)-pow(z,-4)-.4*pow(rc,-10)+pow(rc,-4));
	return .5*(z*z-rmin*rmin)+.4*pow(rmin,-10)-pow(rmin,-4)-.4*pow(rc,-10)+pow(rc,-4);
}


void getUfilterFFT()
{
	for(int i=1;i<=NiLJ;i++)
	{
		double z=dx*i;
		Ufilter[i]=getU(z)*2.*PI;
		Ufilter[N-i]=getU(z)*2.*PI;
	}
	Ufilter[0]=getU(0.)*2.*PI;
	
	
	
	fftw_plan fU = fftw_plan_dft_r2c_1d(N, Ufilter, UfilterFFT, FFTW_ESTIMATE);
	fftw_execute(fU);
	
	
	for(int i=0;i<K;i++)
	{
		UfilterFFT[i][0]*=dx;
		UfilterFFT[i][1]*=dx;
	}
	
  fftw_destroy_plan(fU);
}


void get_mu_LJ()
{
  for(int i=0;i<Nspecies;i++)
  for(int j=0;j<Nspecies;j++)
  {
    mu[i]+= -alphaLJ*eps[i*Nspecies+j]*rhob[j];
  }
  printf("%f\n",mu[0]);
}

void getc1_LJ()
{

  for(int i=0;i<Nspecies;i++)
  for(int j=0;j<Nspecies;j++)
  getc1_LJ_i(i,j,eps[i*Nspecies+j]);

	
}

void getc1_LJ_i(int j1 ,int k1 ,double factor)
{
	for(int i=0;i<N-NiLJ+1;i++)
		rhocopy[i]=rho[IDX(k1,i)]*factor;
			
	for(int i=N-NiLJ+1;i<N;i++)
		rhocopy[i]=0.;
	
	conv_FFT(rhocopy,UfilterFFT,c1_temp);
	add_c(&(c1[IDX(j1,0)]),c1_temp);
	
}
