#include "spline.h"
#include <stdlib.h>
#include <stdio.h>
void
create_spline(double *y,
	      int n,
	      double yp1,
	      double ypn,
	      double *y2)
{
  int i,k;
  double p,qn,sig,un,*u;

  u=(double *)malloc(sizeof(double)*n);
  if(yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]  = (3.0/(1-0))*((y[1]-y[0])/(1-0)-yp1);
  }
  for(i=1;i<n-1;i++) {
    sig=(i-(i-1))/((i+1)-(i-1));
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/((i+1)-i) - (y[i]-y[i-1])/(i-(i-1));
    u[i]=(6.0*u[i]/((i+1)-(i-1))-sig*u[i-1])/p;
  }
  if(ypn>0.99e30)
    qn=un=0;
  else {
    qn=0.5;
    un=(3.0/((n-1)-(n-2)))*(ypn-(y[n-1]-y[n-2])/((n-1)-(n-2)));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for(k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}

double
spline_interpolation(double *ya,
		     double *y2a,
		     int n,
		     double x)
{
  int klo,khi,k;
  double h,b,a,y;

  klo=0;
  khi=n-1;
  while((khi-klo)>1) {
    k=(khi+klo) >> 1;
    if(k>x) 
      khi=k;
    else
      klo=k;
  }
  h=khi-klo;
  if(h==0) {
    printf("Fatal error in spline interpolation.\n");
    exit(1);
  }
  a=(khi-x)/h;
  b=(x-klo)/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  
  return y;
}
