#ifndef INIT_H
#define INIT_H




extern int Nspecies;
extern double  eta;
extern double *rho;
extern double *rhob;
extern double *mu;
extern double *rhonew;
extern double *rhocopy;
extern double *Vext;
extern double *ew;
extern double *c1; //one-body effective potential due to interactions
extern double *c1_temp;
extern double *c1_bulk; //one-body effective potential due to interactions at the bulk Vext=0

extern double *rhobL;
extern double *rhobG;

extern int LJ_exists;
extern int ES_exists;
extern int LG_exists;
extern int Is_polar;
extern double h_target;

int file_exists();
void read_params_geometry();
void initialize_dataframes();
void initialize_vars();
void initialize_Vext();
void initialize_rho();
void read_params_system();
void create_outfname();


#endif
