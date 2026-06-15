#include "../include/init.h"
#include<stdio.h>
#include<math.h>
#include <string.h>
#include <stdlib.h>

#include "../include/globals.h"
#include "../include/utils.h"
#include "../include/functional_fmt.h"
#include "../include/functional_LJ.h"
#include "../include/functional_ES.h"
#include "../include/iterator.h"


int Nspecies=0;
double eta=0.0;
double *ew=NULL;
double *rho=NULL;
double *rhonew=NULL;
double *Vext=NULL;
double *c1=NULL;
double *rhocopy=NULL;
double *c1_temp=NULL;

double *c1_bulk=NULL;// equals mu_ex
double *rhob=NULL;
double *mu=NULL;

double *rhobL=NULL;
double *rhobG=NULL;

int LJ_exists=0;//checks for LJ interactions in the system
int ES_exists=0;//checks for ES interactions in the system
int LG_exists=0;//checks for Liquid-Gas co-existence in the system
int Is_polar=0;//checks whether solvent is polar
double h_target=0.;//target adsorption film height (only if Liquid-Gas interface forms) (Gamma=h*(rhobL-rhobG)) see (Archer et al. 2017) https://doi.org/10.1063/1.4974832

void initialize_dataframes()
{
  //allocate memory
  rho = malloc(Nspecies * N * sizeof(double));
  rhonew = malloc(Nspecies * N * sizeof(double));
  Vext = malloc(Nspecies * N * sizeof(double));
  c1 = malloc(Nspecies * N * sizeof(double));
  
  rhocopy = malloc(N * sizeof(double));
  c1_temp = malloc(N * sizeof(double));
  
  c1_bulk=malloc(Nspecies * sizeof(double));
  rhob=malloc(Nspecies * sizeof(double));
  ew=malloc(Nspecies * sizeof(double));
  mu=malloc(Nspecies * sizeof(double));
  
  if(Is_polar)
  {
    P = malloc(N * sizeof(double));
    Pnew = malloc(N * sizeof(double));
    initialize_val(P,N,0.);
    initialize_val(Pnew,N,0.);
  }
  
  //allocate FMT dataframes
  initialize_FMT_omegas();
  initialize_FMT_n();
  initialize_FMT_dphidn();
  
  //allocate utils dataframes
  initialize_utils();
  
  //allocate iterator dataframes
  initialize_iterator_df();
  
  
  //initialize with zeros
  initialize_val(rho,Nspecies * N,0.);
  initialize_val(rhonew,Nspecies * N,0.);
  initialize_val(Vext,Nspecies * N,0.);
  initialize_val(c1,Nspecies * N,0.);
  
  initialize_val(rhocopy,N,0.);
  initialize_val(c1_temp,N,0.);
  
  initialize_val(c1_bulk,Nspecies,0.);
  initialize_val(rhob,Nspecies,0.);
  initialize_val(ew,Nspecies,0.);
  initialize_val(mu,Nspecies,0.);
  
}


void initialize_Vext()
{
  for(int i=NiR;i<N;i++)
  {
    for(int j=0;j<Nspecies;j++)
      {
        Vext[IDX(j,i)]=ew[j]*(2./15*pow(dx*i,-9)-pow(dx*i,-3));
      }
  }
  
  for(int i=0;i<NiR;i++)
  {
    for(int j=0;j<Nspecies;j++)
      {
        Vext[IDX(j,i)]=1000.;
      }
  }
}

void initialize_rho()
{
  for(int i=0;i<NiR;i++)
  {
    for(int j=0;j<Nspecies;j++)
      {
        rho[IDX(j,i)]=0.;
      }
  }
  
  if(LG_exists==0)
  for(int i=NiR;i<N;i++)
  {
    for(int j=0;j<Nspecies;j++)
      {
        rho[IDX(j,i)]=rhob[j]*exp(-Vext[IDX(j,i)]);
        
      }
  }
  else
  {
    int h_target_N = (int)(h_target/dx);
    //left half init wil rhoL_bulk*exp(-Vext)
    for(int i=NiR;i<h_target_N;i++)
    {
      for(int j=0;j<Nspecies;j++)
        {
          rho[IDX(j,i)]=rhobL[j]*exp(-Vext[IDX(j,i)]);
        }
    }
    
    //right half init wil rhoL_bulk*exp(-Vext)
    for(int i=h_target_N;i<N;i++)
    {
      for(int j=0;j<Nspecies;j++)
        {
          rho[IDX(j,i)]=rhobG[j]*exp(-Vext[IDX(j,i)]);
        }
    }
    
  }  
  //printf("%f %f\n",rho[IDX(0,N/2)],rho[IDX(1,N/2)]);
}

void initialize_vars()
{
  //initialize with system param values ew and rhob
  read_params_system();
  
  initialize_Vext();
  initialize_rho();
  
  get_eta();
  
  //initialize mu
  get_mu_HS();
  
  if(LJ_exists) 
  get_mu_LJ();
 
  
  //initialize c1_bulk
  for(int i=0;i<Nspecies;i++)
  c1_bulk[i]=-mu[i]+log(rhob[i]);
}



void read_params_geometry()
{
  FILE *fp=fopen("../data/params_geometry.txt", "r");
	char line[256];
	while (fgets(line, sizeof(line), fp)) 
	{
    if (sscanf(line, "Nspecies=%d", &Nspecies) == 1) continue;
    if (sscanf(line, "N=%d", &N) == 1) continue;
    if (sscanf(line, "R=%lf", &R) == 1) continue;
    if (sscanf(line, "dx=%lf", &dx) == 1) continue;
    if (sscanf(line, "alpha=%lf", &alpha) == 1) continue;
    if (sscanf(line, "Nbatch=%d", &Nbatch) == 1) continue;
    if (sscanf(line, "Is_polar=%d", &Is_polar) == 1) continue;
  }
  fclose(fp);
  
  sigma=2*R;
  NiR=(int)(R/dx);
  K=(N/2+1);
  Lx=dx*N;
  iend=N-3*NiR;
  PI=3.14159265358979;
}

void read_params_system()
{
    FILE *fp=fopen("../data/params_system.txt", "r");
	  char line[256];
    
    
    while (fgets(line, sizeof(line), fp)) 
    {      
      //read eps
      if (strncmp(line, "eps", 3) == 0)
      {
        LJ_exists=1;
        initialize_LJ_df();
        
        char *p = strchr(line, '[');
          if (p)
          {
            p++;

            for (int i = 0; i < Nspecies*Nspecies; i++)
            {
                eps[i] = strtod(p, &p);

                while (*p == ' ' || *p == ',')
                    p++;
            }
          }

          continue;
        
      }
      
      //read h_target
      if (sscanf(line, "h_target=%lf",&h_target) == 1) continue;
      
      //read rhobL and allocates memory for rhobG
      if (strncmp(line, "rhobL", 5) == 0)
      {
        LG_exists=1;
        
        rhobL=malloc(Nspecies * sizeof(double));
        rhobG=malloc(Nspecies * sizeof(double));
        
        char *p = strchr(line, '[');
          if (p)
          {
            p++;

            for (int i = 0; i < Nspecies; i++)
            {
                rhobL[i] = strtod(p, &p);

                while (*p == ' ' || *p == ',')
                    p++;
            }
          }

          continue;
        
      }
      
      //read rhobG
      if (strncmp(line, "rhobG", 5) == 0)
      {
        char *p = strchr(line, '[');
          if (p)
          {
            p++;

            for (int i = 0; i < Nspecies; i++)
            {
                rhobG[i] = strtod(p, &p);
                rhob[i] = rhobG[i]; //The right end of the box is fixed at rhobG (bulk gas density)
                
                while (*p == ' ' || *p == ',')
                    p++;
            }
          }

          continue;
        
      }
      
      
      
      //read q
      if (strncmp(line, "q", 1) == 0)
      {
        ES_exists=1;
        initialize_ES_df();
        
        char *p = strchr(line, '[');
          if (p)
          {
            p++;

            for (int i = 0; i < Nspecies; i++)
            {
                q[i] = strtod(p, &p);

                while (*p == ' ' || *p == ',')
                    p++;
            }
          }

          continue;
        
      }
      
      //read BC string
      if (sscanf(line, "BC=%s", BC) == 1) continue;
      
      //read lambdaB
      if (sscanf(line, "lambdaB=%lf", &lambdaB) == 1) continue;
      
       //read Vq_L
      if (sscanf(line, "Vq_L=%lf", &Vq_L) == 1) continue;
      
      //read Vq_R
      if (sscanf(line, "Vq_R=%lf", &Vq_R) == 1) continue;
      
      //read p (dipole moment)
      if (sscanf(line, "p=%lf", &p) == 1) continue;
      
      //read ew
      if (strncmp(line, "ew", 2) == 0)
        {
          char *p = strchr(line, '[');
          if (p)
          {
            p++;

            for (int i = 0; i < Nspecies; i++)
            {
                ew[i] = strtod(p, &p);

                while (*p == ' ' || *p == ',')
                    p++;
            }
          }

          continue;
        }
      
      //read rhob
      if (strncmp(line, "rhob", 4) == 0)
        {
          char *p = strchr(line, '[');
          if (p)
          {
            p++;

            for (int i = 0; i < Nspecies; i++)
            {
                rhob[i] = strtod(p, &p);

                while (*p == ' ' || *p == ',')
                    p++;
            }
          }

          continue;
        }
    }
   fclose(fp);

}


