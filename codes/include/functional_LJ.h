#ifndef FUNCTIONAL_LJ_H
#define FUNCTIONAL_LJ_H

#include <fftw3.h>

extern double *eps;
extern double *Ufilter;
extern fftw_complex *UfilterFFT;

extern double rc;
extern double rmin;
extern int NiLJ;

extern double alphaLJ;

void initialize_LJ_df();
double getU(double);
void getUfilterFFT();
void get_mu_LJ();
void getc1_LJ();
void getc1_LJ_i(int,int,double);

#endif
