#include "../include/utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<limits.h>
#include <complex.h>
#include<float.h>
#include<string.h>

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

void poisson_1D(double h,
                int Nt,
                double psi_left,
                double psi_right,
                double *phi,
                double *psi,
                const char *bc)
{
    int M = Nt - 2;   // number of interior points

    if (M <= 0) {
        psi[0] = psi_left;
        psi[Nt-1] = psi_right;
        return;
    }

    if (strcmp(bc, "NN") == 0) {
        fprintf(stderr,
                "poisson_1D: NN boundary conditions are singular.\n");
        return;
    }

    double hh = h * h;

    /* tridiagonal coefficients */
    double *a = malloc(M * sizeof(double));   /* lower */
    double *b = malloc(M * sizeof(double));   /* diagonal */
    double *c = malloc(M * sizeof(double));   /* upper */

    double *bprime = malloc(M * sizeof(double));
    double *dprime = malloc(M * sizeof(double));

    if (!a || !b || !c || !bprime || !dprime) {
        fprintf(stderr, "poisson_1D: memory allocation failed\n");

        free(a);
        free(b);
        free(c);
        free(bprime);
        free(dprime);
        return;
    }

    /* initialize to DD matrix */
    for (int i = 0; i < M; i++) {
        a[i] =  1.0 / hh;
        b[i] = -2.0 / hh;
        c[i] =  1.0 / hh;

        dprime[i] = phi[i + 1];
    }

    /* ---------------------- */
    /* DD : Dirichlet-Dirichlet */
    /* ---------------------- */
    if (strcmp(bc, "DD") == 0) {

        psi[0]    = psi_left;
        psi[Nt-1] = psi_right;

        dprime[0]     -= a[0]     * psi_left;
        dprime[M - 1] -= c[M - 1] * psi_right;
    }

    /* ---------------------- */
    /* DN : Dirichlet-Neumann */
    /* psi'(right)=psi_right  */
    /* ---------------------- */
    else if (strcmp(bc, "DN") == 0) {

        psi[0] = psi_left;

        dprime[0] -= a[0] * psi_left;

        /* modify last row */
        b[M - 1] = -1.0 / hh;
        c[M - 1] =  0.0;

        dprime[M - 1] = phi[Nt - 2] - psi_right / h;
    }

    /* ---------------------- */
    /* ND : Neumann-Dirichlet */
    /* psi'(left)=psi_left    */
    /* ---------------------- */
    else if (strcmp(bc, "ND") == 0) {

        psi[Nt - 1] = psi_right;

        /* modify first row */
        a[0] = 0.0;
        b[0] = -1.0 / hh;

        dprime[0] = phi[1] + psi_left / h;

        dprime[M - 1] -= c[M - 1] * psi_right;
    }

    else {
        fprintf(stderr,
                "poisson_1D: unknown BC type \"%s\"\n",
                bc);

        free(a);
        free(b);
        free(c);
        free(bprime);
        free(dprime);
        return;
    }

    /* Thomas algorithm */

    bprime[0] = b[0];

    for (int i = 1; i < M; i++) {

        double m = a[i] / bprime[i - 1];

        bprime[i] = b[i] - m * c[i - 1];
        dprime[i] -= m * dprime[i - 1];
    }

    psi[Nt - 2] = dprime[M - 1] / bprime[M - 1];

    for (int i = M - 2; i >= 0; i--) {
        psi[i + 1] =
            (dprime[i] - c[i] * psi[i + 2]) / bprime[i];
    }

    /* reconstruct Neumann boundary values */
    if (strcmp(bc, "DN") == 0) {
        psi[Nt - 1] = psi[Nt - 2] + h * psi_right;
    }

    if (strcmp(bc, "ND") == 0) {
        psi[0] = psi[1] - h * psi_left;
    }

    free(a);
    free(b);
    free(c);
    free(bprime);
    free(dprime);
}


/*
void poisson_1D(double h, int Nt,
                double psi_left,
                double psi_right,
                double *phi,
                double *psi)
{
    int M = Nt - 2;  // number of interior points

    if (M <= 0) {
        psi[0] = psi_left;
        psi[Nt-1] = psi_right;
        return;
    }

    // Tridiagonal coefficients
    double a =  1.0 / (h*h);
    double b = -2.0 / (h*h);
    double c =  1.0 / (h*h);

    // Boundary conditions
    psi[0]   = psi_left;
    psi[Nt-1] = psi_right;

    // Allocate temporary arrays
    double *bprime = malloc(M * sizeof(double));
    double *dprime = malloc(M * sizeof(double));

    // RHS (with boundary corrections)
    dprime[0] = phi[1] - a * psi_left;

    for (int i = 1; i < M-1; i++)
        dprime[i] = phi[i+1];

    dprime[M-1] = phi[Nt-2] - c * psi_right;

    // Forward sweep (Thomas algorithm)
    bprime[0] = b;

    for (int i = 1; i < M; i++) {
        double m = a / bprime[i-1];
        bprime[i] = b - m * c;
        dprime[i] -= m * dprime[i-1];
    }

    // Back substitution
    psi[Nt-2] = dprime[M-1] / bprime[M-1];

    for (int i = M-2; i >= 0; i--) {
        psi[i+1] = (dprime[i] - c * psi[i+2]) / bprime[i];
    }

    // Free memory
    free(bprime);
    free(dprime);
}
*/

