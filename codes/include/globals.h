#ifndef GLOBALS_H
#define GLOBALS_H

extern int Nbatch;
extern double PI;
extern double alpha;
extern double dx;
extern int N;
extern int K;
extern double Lx;
extern int iend;

extern double R;
extern int NiR; 
extern double sigma;

extern int count_iter;
extern char *outfname;

#define IDX(i,j) ((i)*N + (j))


#endif
