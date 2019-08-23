#include <math.h>

#include "geometry.h"

#define XX  0
#define YY  1
#define ZZ  2
#define DIM 3

void 
copy_rvec(double * a,double * b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}


void rvec_add(double * a,double * b,double * c)
{
  double x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

void
rvec_sub(double * a,double * b,double * c)
{
  double x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

double norm2(double * a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

void svmul(double a,double *  v1,double *  v2)
{
  v2[XX]=a*v1[XX];
  v2[YY]=a*v1[YY];
  v2[ZZ]=a*v1[ZZ];
}

void unitv(double *  src,double *  dest)
{
  double linv;
  
  linv=1.0/sqrt(norm2(src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}


void oprod(double *  a,double *  b,double *  c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

double iprod(double *  a,double *  b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}



double det(double a[3][3])
{
  return ( a[XX][XX]*(a[YY][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[YY][ZZ])
          -a[YY][XX]*(a[XX][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[XX][ZZ])
          +a[ZZ][XX]*(a[XX][YY]*a[YY][ZZ]-a[YY][YY]*a[XX][ZZ]));
}

void mvmul(double a[3][3],double * src,double * dest)
{
  dest[XX]=a[XX][XX]*src[XX]+a[XX][YY]*src[YY]+a[XX][ZZ]*src[ZZ];
  dest[YY]=a[YY][XX]*src[XX]+a[YY][YY]*src[YY]+a[YY][ZZ]*src[ZZ];
  dest[ZZ]=a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
}

double 
cos_angle_no_table(double * a,double * b)
{
  /* This version does not need the invsqrt lookup table */
  double   cos;
  int      m;
  double   aa,bb,ip,ipa,ipb; /* For accuracy these must be double! */
  
  ip=ipa=ipb=0;
  for(m=0; (m<DIM); m++) {              /* 18           */
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);                 /* 12           */
                                        /* 30 TOTAL     */
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}


void mmul(double a[3][3],double b[3][3],double dest[3][3])
{
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[XX][YY]*b[YY][XX]+a[XX][ZZ]*b[ZZ][XX];
  dest[YY][XX]=a[YY][XX]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[YY][ZZ]*b[ZZ][XX];
  dest[ZZ][XX]=a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
  dest[XX][YY]=a[XX][XX]*b[XX][YY]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[ZZ][YY];
  dest[YY][YY]=a[YY][XX]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[ZZ][YY];
  dest[ZZ][YY]=a[ZZ][XX]*b[XX][YY]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
  dest[XX][ZZ]=a[XX][XX]*b[XX][ZZ]+a[XX][YY]*b[YY][ZZ]+a[XX][ZZ]*b[ZZ][ZZ];
  dest[YY][ZZ]=a[YY][XX]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[YY][ZZ]*b[ZZ][ZZ];
  dest[ZZ][ZZ]=a[ZZ][XX]*b[XX][ZZ]+a[ZZ][YY]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

void m_inv(double src[3][3],double dest[3][3])
{
  double  deter,c,fc;
  
  deter = det(src);
  c     = 1.0/deter;
  fc    = fabs(c);
  
  
  dest[XX][XX]= c*(src[YY][YY]*src[ZZ][ZZ]-src[ZZ][YY]*src[YY][ZZ]);
  dest[XX][YY]=-c*(src[XX][YY]*src[ZZ][ZZ]-src[ZZ][YY]*src[XX][ZZ]);
  dest[XX][ZZ]= c*(src[XX][YY]*src[YY][ZZ]-src[YY][YY]*src[XX][ZZ]);
  dest[YY][XX]=-c*(src[YY][XX]*src[ZZ][ZZ]-src[ZZ][XX]*src[YY][ZZ]);
  dest[YY][YY]= c*(src[XX][XX]*src[ZZ][ZZ]-src[ZZ][XX]*src[XX][ZZ]);
  dest[YY][ZZ]=-c*(src[XX][XX]*src[YY][ZZ]-src[YY][XX]*src[XX][ZZ]);
  dest[ZZ][XX]= c*(src[YY][XX]*src[ZZ][YY]-src[ZZ][XX]*src[YY][YY]);
  dest[ZZ][YY]=-c*(src[XX][XX]*src[ZZ][YY]-src[ZZ][XX]*src[XX][YY]);
  dest[ZZ][ZZ]= c*(src[XX][XX]*src[YY][YY]-src[YY][XX]*src[XX][YY]);
}


double 
calc_dist(double *x, int a1, int a2) 
{
  double dx,dy,dz;
  double r;

  dx=x[3*a1  ]-x[3*a2];
  dy=x[3*a1+1]-x[3*a2+1];
  dz=x[3*a1+2]-x[3*a2+2];

  r=sqrt(dx*dx+dy*dy+dz*dz);

  return r;
}

double 
calc_angle(double *x, int a1, int a2, int a3) 
{
  double theta,ctheta;
  double v1[3],v2[3];

    rvec_sub(x+3*a1,x+3*a2,v1);
    rvec_sub(x+3*a3,x+3*a2,v2);

    ctheta= cos_angle_no_table(v1,v2);
    if(ctheta<=-1.0)
      theta = M_PI;
    else if(ctheta>=1.0)
      theta = 0;
    else
      theta = acos(ctheta);

    return theta;
}

double 
calc_torsion(double *x, int a1, int a2, int a3, int a4) 
{
  double tmp_u[3],tmp_v[3],tmp_w[3],v1[3],tmp_a[3];
  double phi_w, phi,scratch;

    rvec_sub(x+3*a3,x+3*a2,tmp_u);
    unitv(tmp_u,tmp_u);

    rvec_sub(x+3*a1,x+3*a2,tmp_v);
    scratch=iprod(tmp_u,tmp_v);
    svmul(scratch,tmp_u,v1);
    rvec_sub(tmp_v,v1,tmp_v);
    unitv(tmp_v,tmp_v);

    rvec_sub(x+3*a4,x+3*a3,tmp_a);
    scratch=iprod(tmp_a,tmp_u);
    svmul(scratch,tmp_u,v1);
    rvec_sub(tmp_a,v1,tmp_a);
    unitv(tmp_a,tmp_a);
    
    oprod(tmp_v,tmp_u,tmp_w);
    unitv(tmp_w,tmp_w);

    scratch=iprod(tmp_a,tmp_v);

    if(scratch>1.0)
      phi = 0.0;
    else if(scratch<=-1.0)
        phi = M_PI;
    else
        phi = acos(scratch);

    scratch=iprod(tmp_a,tmp_w);

    if(scratch>=1.0)
        phi_w = 0.0;
    else if(scratch<=-1.0)
        phi_w = M_PI;
    else
        phi_w=acos(scratch);

    if(phi_w < M_PI/2.0)
        phi = -phi;
  
    return phi;
}


void rotmatrix_from_vector(double * v, double theta, double rotmat[3][3])
{
  double vu[3];
  double x,y,z;
  unitv(v,vu);

  x=vu[0];
  y=vu[1];
  z=vu[2];

  rotmat[0][0]=x*x+(1-x*x)*cos(theta);
  rotmat[0][1]=x*y*(1-cos(theta))-z*sin(theta);
  rotmat[0][2]=x*z*(1-cos(theta))+y*sin(theta);
  rotmat[1][0]=x*y*(1-cos(theta))+z*sin(theta);
  rotmat[1][1]=y*y+(1-y*y)*cos(theta);
  rotmat[1][2]=y*z*(1-cos(theta))-x*sin(theta);
  rotmat[2][0]=x*z*(1-cos(theta))-y*sin(theta);
  rotmat[2][1]=y*z*(1-cos(theta))+x*sin(theta);
  rotmat[2][2]=z*z+(1-z*z)*cos(theta);

}

void rotvector_from_matrix(double rotmat[3][3], double * v, double *theta)
{
  double costheta,sgn;

  costheta = 0.5*(rotmat[0][0]+rotmat[1][1]+rotmat[2][2]-1);

  if(costheta<=-1.0)
    *theta=M_PI;
  else if(costheta>=1.0)
    *theta=0;
  else
    *theta = acos( costheta );

  if(sin(*theta)>0)
    sgn = 1;
  else
    sgn = -1;
  
  if(*theta<1e-7) {
    v[0]=v[1]=v[2]=1.0;
    *theta=0;
  } else {
    v[0] = ( rotmat[2][1]-rotmat[1][2] ) * sgn;
    v[1] = ( rotmat[0][2]-rotmat[2][0] ) * sgn;
    v[2] = ( rotmat[1][0]-rotmat[0][1] ) * sgn;
  }
  unitv(v,v);

}

void calc_cog(double * x1, double * x2, double * x3, double * cog)
{
  int d;
  
  cog[XX]=cog[YY]=cog[ZZ]=0;
  
  for(d=0;d<DIM;d++) 
    cog[d] += (x1[d] + x2[d] + x3[d])/3.0;
}


void 
calc_fit3(double * a1, double * a2, double * a3, 
	 double * b1, double * b2, double * b3,
	 double * cogA, double * cogB, 
	 double * rotvec, double *rotangle) 
{
  double rotmat[3][3];
  double tm1[3][3];
  double tm2[3][3];
  double vA[3][3],vB[3][3];
  double tv[3];
  double scratch;
  
  
  calc_cog(a1,a2,a3,cogA);
  calc_cog(b1,b2,b3,cogB);

  rvec_sub(a2,a1,vA[0]);
  unitv(vA[0],vA[0]);
  rvec_sub(a3,a1,vA[1]);
  scratch=iprod(vA[0],vA[1]);
  svmul(scratch,vA[0],tv);
  rvec_sub(vA[1],tv,vA[1]);
  unitv(vA[1],vA[1]);

  oprod(vA[0],vA[1],vA[2]);

  rvec_sub(b2,b1,vB[0]);
  unitv(vB[0],vB[0]);
  rvec_sub(b3,b1,vB[1]);
  scratch=iprod(vB[0],vB[1]);
  svmul(scratch,vB[0],tv);
  rvec_sub(vB[1],tv,vB[1]);
  unitv(vB[1],vB[1]);

  oprod(vB[0],vB[1],vB[2]);

  m_inv(vA,tm1);  
  mmul(tm1,vB,tm2);
  m_inv(tm2,rotmat);
  rotvector_from_matrix(rotmat,rotvec,rotangle);

}
