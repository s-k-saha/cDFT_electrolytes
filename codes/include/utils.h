#ifndef UTILS_H
#define UTILS_H

#include <fftw3.h>

extern fftw_complex *fft_f;
extern fftw_complex *fft_h;

void initialize_val(double*,int,double);
void initialize_utils();
void conv_FFT(double*, fftw_complex*,double*);
void add_c(double*, double*);
double getMax(double*,int);
void poisson_1D(double,int,double,double,double*,double*,const char*);
void psi_calculator(double,int,double,double,double*,double*);

#endif
