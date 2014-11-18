/**
 ** routines for allocating, freeing and copying 2-D arrays
 **/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "arrays.h"

void	allocShortArray(short ***f, int N1, int N2)
{
	int		i;
	
	(*f) = (short **) calloc(N2, sizeof(short *));
	if (!*f)
		error_exit("allocShortArray: memory allocation failure\n");
	for (i = 0; i < N2; i++)
		{
			(*f)[i] = (short *) calloc(N1, sizeof(short));
			if (!(*f)[i])
				error_exit("allocShortArray: memory allocation failure\n");
		}
}

void	allocIntArray(int ***f, int N1, int N2)
{
	int		i;
	
	(*f) = (int **) calloc(N2, sizeof(int *));
	if (!*f)
	  error_exit("allocIntArray: memory allocation failure\n");
	for (i = 0; i < N2; i++)
		{
		  (*f)[i] = (int *) calloc(N1, sizeof(int));
		  if (!(*f)[i])
		    error_exit("allocIntArray: memory allocation failure\n");
		}
}


void	allocFloatArray(float ***f, int N1, int N2)
{
  int		i;
  
  (*f) = (float **) calloc(N2, sizeof(float *));
  if (!*f)
    error_exit("allocFloatArray: memory allocation failure\n");
  for (i = 0; i < N2; i++)
    {
      (*f)[i] = (float *) calloc(N1, sizeof(float));
      if (!(*f)[i])
	error_exit("allocFloatArray: memory allocation failure\n");
    }
}



void	allocDoubleArray(double ***f, int N1, int N2)
{
	int		i;
	
	(*f) = (double **) calloc(N2, sizeof(double *));
	if (!*f)
		error_exit("allocDoubleArray: memory allocation failure\n");
	for (i = 0; i < N2; i++)
		{
			(*f)[i] = (double *) calloc(N1, sizeof(double));
			if (!(*f)[i])
				error_exit("allocDoubleArray: memory allocation failure\n");
		}
}





void	freeShortArray(short **f, int N1, int N2)
{
	int	i;
	
	for (i = 0; i < N2; i++)
		free(f[i]);
	free(f);
}


void	freeIntArray(int **f, int N1, int N2)
{
	int	i;
	
	for (i = 0; i < N2; i++)
	  free(f[i]);
	free(f);
}





void	freeFloatArray(float **f, int N1, int N2)
{
	int	i;
	
	for (i = 0; i < N2; i++)
		free(f[i]);
	free(f);
}


void	freeDoubleArray(double **f, int N1, int N2)
{
	int	i;
	
	for (i = 0; i < N2; i++)
		free(f[i]);
	free(f);
}





void	copyFloatToShort(float **fsrc, short **fdst, int N1, int N2)
{
	int	i, j;
	
	for (i = 0; i < N2; i++)
		for (j = 0; j < N1; j++)
			fdst[i][j] = fsrc[i][j];
}





void	copyShortToFloat(short **fsrc, float **fdst, int N1, int N2)
{
	int	i, j;
	
	for (i = 0; i < N2; i++)
		for (j = 0; j < N1; j++)
			fdst[i][j] = fsrc[i][j];
}



