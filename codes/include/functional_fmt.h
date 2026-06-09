#ifndef FUNCTIONAL_FMT_H
#define FUNCTIONAL_FMT_H

#include <fftw3.h>

extern fftw_complex  *omega3;
extern fftw_complex  *omega2;
extern fftw_complex  *omega1;
extern fftw_complex  *omega0;

extern fftw_complex *omega1v;
extern fftw_complex *omega2v;

extern double *n0;
extern double *n1;
extern double *n2;
extern double *n3;
extern double *n1v;
extern double *n2v;

extern double *dphidn0;
extern double *dphidn1;
extern double *dphidn2;
extern double *dphidn3;
extern double *dphidn1v;
extern double *dphidn2v;


void get_eta();
void get_mu_HS();
void initialize_FMT_omegas();
void initialize_FMT_n();
void initialize_FMT_dphidn();

void getn();//FMT weighted densities
void getc1_fmt();




#endif
