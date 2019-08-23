#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rmsd.h"


double
calcrmsd(double *     x1,
     double *     x2,
     int          natoms)
{
  int i;
  double dx,dy,dz;
  double rms;

  rms = 0;
  for(i=0;i<natoms;i++) {
    dx   = x1[3*i]   - x2[3*i];
    dy   = x1[3*i+1] - x2[3*i+1];
    dz   = x1[3*i+2] - x2[3*i+2];
    rms += dx*dx + dy*dy + dz*dz;
  }
  return sqrt(rms/natoms);
}



double
calcrmsdindex(double *     x1,
	      double *     x2,
	      int *        index,
	      int          nindex)
{
  int i;
  double dx,dy,dz;
  double rms;

  rms = 0;
  for(i=0;i<nindex;i++) {
    dx   = x1[3*index[i]]   - x2[3*index[i]];
    dy   = x1[3*index[i]+1] - x2[3*index[i]+1];
    dz   = x1[3*index[i]+2] - x2[3*index[i]+2];
    rms += dx*dx + dy*dy + dz*dz;
  }
  return sqrt(rms/nindex);
}




#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);
	
static void 
jacobi(double **a,int n,double d[],double **v,int *nrot)
{
  int j,i;
  int iq,ip;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b=(double *)malloc(sizeof(double)*n);
  z=(double *)malloc(sizeof(double)*n);

  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
 }
  for (ip=0; ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1; i<=500; i++) {
    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      free(z);
      free(b);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
            && fabs(d[iq])+g == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if (fabs(h)+g == fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0; j<ip; j++) {
            ROTATE(a,j,ip,j,iq)
	  }
          for (j=ip+1; j<iq; j++) {
            ROTATE(a,ip,j,j,iq)
            }
          for (j=iq+1; j<n; j++) {
            ROTATE(a,ip,j,iq,j)
            }
          for (j=0; j<n; j++) {
            ROTATE(v,j,ip,j,iq)
            }
          ++(*nrot);
        }
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] +=  z[ip];
      d[ip]  =  b[ip];
      z[ip]  =  0.0;
    }
  }
  fprintf(stderr,"Error: Too many iterations in routine JACOBI\n");
}


static void 
calc_fit_R(int natoms,double *xp,double *x,double R[3][3])
{
  int    c,r,n,j,i,irot;
  double **omega,**om;
  double d[2*3],xnr,xpc;
  double vh[3][3];
  double vk[3][3];
  double u[3][3];

  int    index;
  double   max_d;

  omega=(double **)malloc(sizeof(double *)*2*3);
  om=(double **)malloc(sizeof(double *)*2*3);
  
  for(i=0; i<2*3; i++) {
    omega[i]=(double *)malloc(sizeof(double)*2*3);
    om[i]=(double *)malloc(sizeof(double)*2*3);
  }
  
  for(i=0; i<2*3; i++) {
    d[i]=0;
    for(j=0; j<2*3; j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      u[i][j]=0;

  for(n=0;(n<natoms);n++)
      for(c=0; (c<3); c++) {
	xpc=xp[n*3+c];
	for(r=0; (r<3); r++) {
	  xnr=x[n*3+r];
	  u[c][r]+=xnr*xpc;
	}
      }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<2*3; r++)
    for(c=0; c<=r; c++)
      if (r>=3 && c<3) {
        omega[r][c]=u[r-3][c];
        omega[c][r]=u[r-3][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }
  
  /*determine h and k*/
  jacobi(omega,2*3,d,om,&irot);
  
  index=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0; j<2; j++) {
    max_d=-1000;
    for(i=0; i<2*3; i++)
      if (d[i]>max_d) {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0; i<3; i++) {
      vh[j][i]=M_SQRT2*om[i][index];
      vk[j][i]=M_SQRT2*om[i+3][index];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  

  vh[2][0]=vh[0][1]*vh[1][2]-vh[0][2]*vh[1][1];
  vh[2][1]=vh[0][2]*vh[1][0]-vh[0][0]*vh[1][2];
  vh[2][2]=vh[0][0]*vh[1][1]-vh[0][1]*vh[1][0];

  vk[2][0]=vk[0][1]*vk[1][2]-vk[0][2]*vk[1][1];
  vk[2][1]=vk[0][2]*vk[1][0]-vk[0][0]*vk[1][2];
  vk[2][2]=vk[0][0]*vk[1][1]-vk[0][1]*vk[1][0];


  /*determine R*/
  for(r=0; r<3; r++)
    for(c=0; c<3; c++)
      R[r][c] = vk[0][r]*vh[0][c] +
	        vk[1][r]*vh[1][c] +
	        vk[2][r]*vh[2][c];

  for(i=0; i<2*3; i++) {
    free(omega[i]);
    free(om[i]);
  }
  free(omega);
  free(om);
}


static void 
calc_fit_R_index(int *index,int nindex,double *xp,double *x,double R[3][3])
{
  int    c,r,n,j,i,irot;
  double **omega,**om;
  double d[2*3],xnr,xpc;
  double vh[3][3];
  double vk[3][3];
  double u[3][3];

  int    idx;
  double   max_d;

  omega=(double **)malloc(sizeof(double *)*2*3);
  om=(double **)malloc(sizeof(double *)*2*3);
  
  for(i=0; i<2*3; i++) {
    omega[i]=(double *)malloc(sizeof(double)*2*3);
    om[i]=(double *)malloc(sizeof(double)*2*3);
  }
  
  for(i=0; i<2*3; i++) {
    d[i]=0;
    for(j=0; j<2*3; j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      u[i][j]=0;

  for(n=0;(n<nindex);n++)
      for(c=0; (c<3); c++) {
	xpc=xp[n*3+c];
	for(r=0; (r<3); r++) {
	  xnr=x[index[n]*3+r];
	  u[c][r]+=xnr*xpc;
	}
      }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<2*3; r++)
    for(c=0; c<=r; c++)
      if (r>=3 && c<3) {
        omega[r][c]=u[r-3][c];
        omega[c][r]=u[r-3][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }
  
  /*determine h and k*/
  jacobi(omega,2*3,d,om,&irot);
  
  idx=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0; j<2; j++) {
    max_d=-1000;
    for(i=0; i<2*3; i++)
      if (d[i]>max_d) {
        max_d=d[i];
        idx=i;
      }
    d[idx]=-10000;
    for(i=0; i<3; i++) {
      vh[j][i]=M_SQRT2*om[i][idx];
      vk[j][i]=M_SQRT2*om[i+3][idx];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  

  vh[2][0]=vh[0][1]*vh[1][2]-vh[0][2]*vh[1][1];
  vh[2][1]=vh[0][2]*vh[1][0]-vh[0][0]*vh[1][2];
  vh[2][2]=vh[0][0]*vh[1][1]-vh[0][1]*vh[1][0];

  vk[2][0]=vk[0][1]*vk[1][2]-vk[0][2]*vk[1][1];
  vk[2][1]=vk[0][2]*vk[1][0]-vk[0][0]*vk[1][2];
  vk[2][2]=vk[0][0]*vk[1][1]-vk[0][1]*vk[1][0];


  /*determine R*/
  for(r=0; r<3; r++)
    for(c=0; c<3; c++)
      R[r][c] = vk[0][r]*vh[0][c] +
	        vk[1][r]*vh[1][c] +
	        vk[2][r]*vh[2][c];

  for(i=0; i<2*3; i++) {
    free(omega[i]);
    free(om[i]);
  }
  free(omega);
  free(om);
}


void
fitrmsd(double *  ref_x,
	double *  x,
	int       natoms)
{
  double *tmprefx;
  double rcx,rcy,rcz;
  double cx,cy,cz;
  double R[3][3];
  double xo,yo,zo;
  int i,j;

  tmprefx = (double *)malloc(sizeof(double)*3*natoms);
  
  rcx=rcy=rcz=0;
  for(i=0;i<natoms;i++) {
    rcx += ref_x[3*i];
    rcy += ref_x[3*i+1];
    rcz += ref_x[3*i+2];
  }
  rcx/=natoms;
  rcy/=natoms;
  rcz/=natoms;

  for(i=0;i<natoms;i++) {
    tmprefx[3*i]   = ref_x[3*i]   - rcx;
    tmprefx[3*i+1] = ref_x[3*i+1] - rcy;
    tmprefx[3*i+2] = ref_x[3*i+2] - rcz;
  }

  cx=cy=cz=0;
  for(i=0;i<natoms;i++) {
    cx += x[3*i];
    cy += x[3*i+1];
    cz += x[3*i+2];
  }
  cx/=natoms;
  cy/=natoms;
  cz/=natoms;
  for(i=0;i<natoms;i++) {
    x[3*i]   -= cx;
    x[3*i+1] -= cy;
    x[3*i+2] -= cz;
  }

  /* Both centered in origin now */
  
  /* Calc fit */
  calc_fit_R(natoms,tmprefx,x,R);

  /* Rotate and translate */
  for(j=0; j<natoms; j++) {
    xo = x[3*j];
    yo = x[3*j+1];
    zo = x[3*j+2];
    x[3*j]   = R[0][0]*xo + R[0][1]*yo + R[0][2]*zo + rcx;
    x[3*j+1] = R[1][0]*xo + R[1][1]*yo + R[1][2]*zo + rcy;
    x[3*j+2] = R[2][0]*xo + R[2][1]*yo + R[2][2]*zo + rcz;
  }  


}
  
  

void
fitrmsdindex(double *  ref_x,
	     double *  x,
	     int *     index,
	     int       nindex)
{
  double *tmprefx;
  double rcx,rcy,rcz;
  double cx,cy,cz;
  double R[3][3];
  double xo,yo,zo;
  int i,j;

  tmprefx = (double *)malloc(sizeof(double)*3*nindex);
  
  rcx=rcy=rcz=0;
  for(i=0;i<nindex;i++) {
    rcx += ref_x[3*index[i]];
    rcy += ref_x[3*index[i]+1];
    rcz += ref_x[3*index[i]+2];
  }
  rcx/=nindex;
  rcy/=nindex;
  rcz/=nindex;

  for(i=0;i<nindex;i++) {
    tmprefx[3*i]   = ref_x[3*index[i]]   - rcx;
    tmprefx[3*i+1] = ref_x[3*index[i]+1] - rcy;
    tmprefx[3*i+2] = ref_x[3*index[i]+2] - rcz;
  }

  cx=cy=cz=0;
  for(i=0;i<nindex;i++) {
    cx += x[3*index[i]];
    cy += x[3*index[i]+1];
    cz += x[3*index[i]+2];
  }
  cx/=nindex;
  cy/=nindex;
  cz/=nindex;
  for(i=0;i<nindex;i++) {
    x[3*index[i]]   -= cx;
    x[3*index[i]+1] -= cy;
    x[3*index[i]+2] -= cz;
  }

  /* Both centered in origin now */
  
  calc_fit_R_index(index,nindex,tmprefx,x,R);

  /* Rotate and translate */
  for(j=0; j<nindex; j++) {
    xo = x[3*index[j]];
    yo = x[3*index[j]+1];
    zo = x[3*index[j]+2];
    x[3*index[j]]   = R[0][0]*xo + R[0][1]*yo + R[0][2]*zo + rcx;
    x[3*index[j]+1] = R[1][0]*xo + R[1][1]*yo + R[1][2]*zo + rcy;
    x[3*index[j]+2] = R[2][0]*xo + R[2][1]*yo + R[2][2]*zo + rcz;
  }  

}
  
  
