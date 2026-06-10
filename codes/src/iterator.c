#include "../include/iterator.h"

#include<stdio.h>
#include<math.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "../include/globals.h"
#include "../include/init.h"
#include "../include/utils.h"
#include "../include/functional_fmt.h"
#include "../include/functional_ES.h"
#include "../include/functional_LJ.h"


double *dev=NULL;

void initialize_iterator_df()
{
  dev=malloc(Nspecies*sizeof(double));
  
  for(int i=0;i<Nspecies;i++)
  dev[i]=0.;
}

void rhocpy()
{
  for(int i=0;i<N;i++)
  rhocopy[i]=0.;
  
  for(int j=0;j<Nspecies;j++)
  {
    for(int i=0;i<N-NiR+1;i++)
		  rhocopy[i]+=rho[IDX(j,i)];
	}
	
	for(int i=N-NiR+1;i<N;i++)
		rhocopy[i]=0.0;
}

void iterate()
{
  initialize_val(c1,Nspecies*N,0.0);//clean c1 before each iteration
  rhocpy();
  getn();
  
  getc1_fmt();
  if(LJ_exists)
  getc1_LJ();
  
  if(ES_exists)
  getc1_ES();
  
  if(LG_exists==0)
  for(int j=0;j<Nspecies;j++)
  for(int i=0;i<iend;i++)
	{	
		rhonew[IDX(j,i)]=rhob[j]*exp(-Vext[IDX(j,i)]+c1[IDX(j,i)]-c1_bulk[j]);
	}
	else
	for(int j=0;j<Nspecies;j++)
  for(int i=0;i<iend;i++)
	{	
		rhonew[IDX(j,i)]=exp(-Vext[IDX(j,i)]+c1[IDX(j,i)]+mu[j]);
	}
	
	if(LG_exists==0)
	for(int j=0;j<Nspecies;j++)
  {
    dev[j]=0.;
    for(int i=0;i<iend;i++)
	   dev[j]+=fabs(rhonew[IDX(j,i)]-rho[IDX(j,i)]);
	  dev[j]*=dx;
	}
	
	if(LG_exists==0)
  for(int i=0;i<iend;i++)
	  for(int j=0;j<Nspecies;j++)
		  rho[IDX(j,i)]=(1-alpha)*rho[IDX(j,i)]+alpha*rhonew[IDX(j,i)];
  else
  {
    for(int j=0;j<Nspecies;j++)
    {
      dev[j]=0.;
      double Gamma_new=0.;
      for(int i=NiR;i<iend;i++)
      {
		    rhonew[IDX(j,i)]=(1-alpha)*rho[IDX(j,i)]+alpha*rhonew[IDX(j,i)];
		    Gamma_new+=(rhonew[IDX(j,i)]-rhobG[j]);
      }
      
      Gamma_new*=dx;
      double factor=(h_target*(rhobL[j]-rhobG[j]))/Gamma_new;
        
      for(int i=NiR;i<iend;i++)
      {
		    rhonew[IDX(j,i)]=factor*(rhonew[IDX(j,i)]-rhobG[j])+ rhobG[j];
		    dev[j]+=fabs(rhonew[IDX(j,i)]-rho[IDX(j,i)]);
		    rho[IDX(j,i)]=rhonew[IDX(j,i)];
      }
      dev[j]*=dx;
    }
  }
    
  
}


void write_rho(double elapsed, int iter)
{
    char fname[1024];
    char outname[1024];
    char tmpfile[1100];
    char template_name[250] = "rhob%d_%f";
    char tempname[256];

    /* Construct output filename */
    sprintf(fname, "../data/rho1Ddx%fL%f", dx, Lx);

    for (int i = 0; i < Nspecies; i++)
    {
        sprintf(tempname, template_name, i + 1, rhob[i]);
        strcat(fname, tempname);
    }

    strcpy(template_name, "ew%d_%f");
    for (int i = 0; i < Nspecies; i++)
    {
        sprintf(tempname, template_name, i + 1, ew[i]);
        strcat(fname, tempname);
    }

    if (ES_exists)
    {
        strcpy(template_name, "lambdaB%fVq%fBC%s");
        sprintf(tempname, template_name, lambdaB, Vq, BC);
        strcat(fname, tempname);
    }
    
    if (LG_exists)
    {
        strcpy(template_name, "h_target%f");
        sprintf(tempname, template_name, h_target);
        strcat(fname, tempname);
    }
    
    if (LJ_exists)
        strcat(fname, "_LJ");

    if (ES_exists)
        strcat(fname, "_ES");
    
    if (LG_exists)
        strcat(fname, "_LG");
        
    strcat(fname, ".dat");
    
    strcpy(outname,fname);
    /* Temporary file name */
    snprintf(tmpfile, sizeof(tmpfile), "%s.tmp", fname);

    /* Open temporary file */
    FILE *F = fopen(tmpfile, "w");
    if (F == NULL)
    {
        perror("fopen");
        return;
    }

    /* Write data */
    for (int i = 0; i < N; i++)
    {
        fprintf(F, "%d %f ", i, dx * i);

        for (int j = 0; j < Nspecies; j++)
            fprintf(F, "%f %f ",
                    rho[IDX(j, i)],
                    c1[IDX(j, i)]);

        if (ES_exists)
            fprintf(F, "%f ", ((i >= NiR) && (i <= iend) ) ? psi[i-NiR] : 0.0);

        fprintf(F, "\n");
    }

    fprintf(F,
            "------------------\n%d x %d cycles time: %f s\n",
            Nbatch, iter, elapsed);

    fname[0] = '\0';
    strcpy(template_name, "mu[%d] :%f; ");

    for (int i = 0; i < Nspecies; i++)
    {
        sprintf(tempname, template_name, i, mu[i]);
        strcat(fname, tempname);
    }

    fprintf(F, "%s\n", fname);
    
    
    if (LG_exists)
    {
      fname[0] = '\0';
      strcpy(template_name, "h_target :%f; ");
      sprintf(tempname, template_name, h_target);
      strcat(fname, tempname);
      fprintf(F, "%s\n", fname);
    }
    
    
    if (LJ_exists)
    {
        fname[0] = '\0';
        strcpy(template_name, "eps[%d][%d] :%f; ");

        for (int i = 0; i < Nspecies; i++)
        {
            for (int j = i; j < Nspecies; j++)
            {
                sprintf(tempname,
                        template_name,
                        i,
                        j,
                        eps[i * Nspecies + j]);
                strcat(fname, tempname);
            }
        }

        fprintf(F, "%s\n", fname);
    }

    fname[0] = '\0';
    strcpy(template_name, "dev[%d] :%f; ");

    for (int i = 0; i < Nspecies; i++)
    {
        sprintf(tempname, template_name, i, dev[i]);
        strcat(fname, tempname);
    }

    fprintf(F, "%s\n", fname);

    /* Make sure everything reaches the OS buffers */
    fflush(F);

    if (fclose(F) != 0)
    {
        perror("fclose");
        remove(tmpfile);
        return;
    }

    /*
     * Atomically replace old file with new file.
     * If the program dies before this point,
     * the original file is still intact.
     */
    if (rename(tmpfile, outname) != 0)
    {
        perror("rename");
        remove(tmpfile);
        return;
    }
}



