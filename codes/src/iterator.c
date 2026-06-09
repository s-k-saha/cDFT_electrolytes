#include "../include/iterator.h"

#include<stdio.h>
#include<math.h>
#include <string.h>
#include <stdlib.h>

#include "../include/globals.h"
#include "../include/init.h"
#include "../include/utils.h"
#include "../include/functional_fmt.h"
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
  
  for(int j=0;j<Nspecies;j++)
  for(int i=0;i<iend;i++)
	{	
		rhonew[IDX(j,i)]=rhob[j]*exp(-Vext[IDX(j,i)]+c1[IDX(j,i)]-c1_bulk[j]);
	}
	
	for(int j=0;j<Nspecies;j++)
  {
    dev[j]=0.;
    for(int i=0;i<iend;i++)
	   dev[j]+=fabs(rhonew[IDX(j,i)]-rho[IDX(j,i)]);
	  dev[j]*=dx;
	}
	
	
  for(int i=0;i<iend;i++)
	  for(int j=0;j<Nspecies;j++)
		  rho[IDX(j,i)]=(1-alpha)*rho[IDX(j,i)]+alpha*rhonew[IDX(j,i)];
}

void write_rho(double elapsed,int iter)
{
  //naming the output datafile
  char fname[1024];
  char template_name[250]="rhob%d_%f";
  char tempname[128];
  
	sprintf(fname,"../data/rho1Dew%fdx%fL%f",ew,dx,Lx);
	for(int i=0;i<Nspecies;i++)
	{
	  sprintf(tempname,template_name,i+1,rhob[i]);
	  strcat(fname,tempname);
	}
	
	if(LJ_exists)
	strcat(fname,"_LJ");
	
	strcat(fname,".dat");
	
	
	//writing to file
	FILE *F=fopen(fname,"w");
	for(int i=0;i<N;i++)
	{
	  fprintf(F,"%d %f ",i,dx*i);
	  for(int j=0;j<Nspecies;j++)
		  fprintf(F,"%f %f ",rho[IDX(j,i)],c1[IDX(j,i)]);
		fprintf(F,"\n");
	}
	
	fprintf(F,"------------------\n%d x %d cycles time: %f s\n",Nbatch,iter,elapsed);
	
	strcpy(fname,"");
	strcpy(template_name,"mu[%d] :%f; ");
	for(int i=0;i<Nspecies;i++)
	{
	  sprintf(tempname,template_name,i,mu[i]);
	  strcat(fname,tempname);
	}
	fprintf(F,fname);
	fprintf(F,"\n");
	
	
	if(LJ_exists)
	{
	  strcpy(fname,"");
	  strcpy(template_name,"eps[%d][%d] :%f; ");
	  for(int i=0;i<Nspecies;i++)
	  for(int j=i;j<Nspecies;j++)
	  {
	    sprintf(tempname,template_name,i,j,eps[i*Nspecies+j]);
	    strcat(fname,tempname);
	  }
	  fprintf(F,fname);
	  fprintf(F,"\n");
	}
	
	strcpy(fname,"");
	strcpy(template_name,"dev[%d] :%f; ");
	for(int i=0;i<Nspecies;i++)
	{
	  sprintf(tempname,template_name,i,dev[i]);
	  strcat(fname,tempname);
	}
	fprintf(F,fname);
	fprintf(F,"\n");
	
	
	
	
	fclose(F);
}

