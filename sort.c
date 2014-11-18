#include <string.h>
#include <stdlib.h>

void hpsort(n,ra,ind)
float *ra;
int *ind;
long n;
{
  long i,ir,j,l,irra;
  float rra;
  int *tempind;
  
  tempind = (int *)calloc(n+1,sizeof(int));
  for (i=1;i<n+1;i++)
    tempind[i] = i-1;

  if (n >= 2) {
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
      if (l > 1) {
	irra=tempind[--l];
	rra=ra[irra];
      } else {
	irra=tempind[ir];
	rra=ra[irra];
	tempind[ir]=tempind[1];
	if (--ir == 1) {
	  tempind[1]=irra;
	  break;
	}
      }
      i=l;
      j=l+l;
      while (j <= ir) {
	if (j < ir && ra[tempind[j]] < ra[tempind[j+1]]) j++;
	if (rra < ra[tempind[j]]) {
	  tempind[i]=tempind[j];
	  i=j;
	  j <<= 1;
	} else j=ir+1;
      }
      tempind[i]=irra;
    }
  }

  for (i=1;i<n+1;i++)
    ind[i-1] = tempind[i];
  free(tempind);
}

extern void hpsortstr(long n, char **sa, int *ind)
{
  long i,ir,j,l,irra;
  char *ssa;
  int *tempind;
  
  tempind = (int *)calloc(n+1,sizeof(int));
  for (i=1;i<n+1;i++)
    tempind[i] = i-1;

  if (n >= 2) {
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
      if (l > 1) {
	irra=tempind[--l];
	ssa=sa[irra];
      } else {
	irra=tempind[ir];
	ssa=sa[irra];
	tempind[ir]=tempind[1];
	if (--ir == 1) {
	  tempind[1]=irra;
	  break;
	}
      }
      i=l;
      j=l+l;
      while (j <= ir) {
	if (j < ir && strcmp(sa[tempind[j]],sa[tempind[j+1]])<0) j++;
	if (strcmp(ssa,sa[tempind[j]])<0) {
	  tempind[i]=tempind[j];
	  i=j;
	  j <<= 1;
	} else j=ir+1;
      }
      tempind[i]=irra;
    }
  }

  for (i=1;i<n+1;i++)
    ind[i-1] = tempind[i];
  free(tempind);
}
