#include <math.h>

extern void chdecomp(float **a, int n, float p[], int *err)
{
  int i,j,k;
  float sum;
  
  *err = 0;
  for (i=0;i<n;i++) {
    for (j=i;j<n;j++) {
      for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
      if (i == j) {
	if (sum <= 0.0) {
	  *err = 1;
	  return;
	}
	p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
}