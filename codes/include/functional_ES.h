#ifndef FUNCTIONAL_ES_H
#define FUNCTIONAL_ES_H


extern double lambdaB;
extern double Vq_L;
extern double Vq_R;
extern double *q;
extern double *phi;
extern double *psi;
extern int N_ES;
extern char *BC; 

extern double *E_exec;
extern double *P;
extern double *Pnew;
extern double p;

void initialize_ES_df();
void getc1_ES();

#endif
