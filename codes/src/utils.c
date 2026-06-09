#include "../include/utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<limits.h>
#include <complex.h>
#include<float.h>

#include "../include/globals.h"
#include "../include/init.h"


fftw_complex *fft_f=NULL;	
fftw_complex *fft_h=NULL;
	
void initialize_val(double *temp,int Nt,double val)
{
  for(int i=0;i<Nt;i++)
  temp[i]=val;
}

void initialize_utils()
{
  fft_f=fftw_alloc_complex((size_t) K);
  fft_h=fftw_alloc_complex((size_t) K);
}

double getMax(double *arr,int Nt)
{
  double maxel=-DBL_MAX;
  for(int i=0;i<Nt;i++)
    maxel=(arr[i]>=maxel)? arr[i]:maxel;
  return maxel;
}

void conv_FFT(double *f, fftw_complex *fft_g,double *h2)
{
	
	fftw_plan ff = fftw_plan_dft_r2c_1d(N, f, fft_f, FFTW_ESTIMATE);
  fftw_execute(ff);
	
  for (int i = 0; i < K; ++i) {
        double a = fft_f[i][0], b = fft_f[i][1];
        double c = fft_g[i][0], d = fft_g[i][1]; 
        fft_h[i][0] = a*c - b*d;
        fft_h[i][1] = a*d + b*c;
    }
  
  
	fftw_plan fh= fftw_plan_dft_c2r_1d(N,fft_h,h2,FFTW_ESTIMATE);
	fftw_execute(fh);
	
	
	fftw_destroy_plan(ff);
	fftw_destroy_plan(fh);
	
	
	for(int i=0;i<N;i++)
	h2[i]=(h2[i]/N);
	
	//fftw_free(fft_f);
  //fftw_free(fft_h);
}

void add_c(double *c_1, double *c_temp)
{
	for(int i=0;i<N;i++)
	c_1[i]-=c_temp[i];
}

