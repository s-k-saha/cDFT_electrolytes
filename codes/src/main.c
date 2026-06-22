#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include<limits.h>

#include "../include/globals.h"
#include "../include/init.h"
#include "../include/utils.h"
#include "../include/functional_fmt.h"
#include "../include/functional_LJ.h"
#include "../include/functional_ES.h"
#include "../include/iterator.h"


int main(int argc,char *argv[])
{
  read_params_geometry();
  initialize_dataframes();
  initialize_vars();
  
  clock_t start = clock();clock_t end;
	double elapsed;
	
	
	//repeated iterations
	for(int i=0;i<INT_MAX;++i)
	{	
		for(int j=0;j<Nbatch;j++)
			iterate();
		end = clock();
		elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		write_rho(elapsed);
		printf("%d x %d batches; time: %f s; max dev: %lf\n",count_iter/Nbatch,Nbatch,elapsed,getMax(dev,Nspecies));
	}

}
