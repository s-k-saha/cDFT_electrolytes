#ifndef FUNCTIONAL_ES_H
#define FUNCTIONAL_ES_H


extern double lambdaB;
extern double Vq;
extern double *q;
extern double *phi;
extern double *psi;
extern int N_ES;
extern char *BC; 

void initialize_ES_df();
void getc1_ES();

#endif
