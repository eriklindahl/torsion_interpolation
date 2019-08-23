#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

double 
calc_dist(double *x, int a1, int a2);


double 
calc_angle(double *x, int a1, int a2, int a3);


double 
calc_torsion(double *x, int a1, int a2, int a3, int a4);

void
rvec_sub(double * a,double * b,double * c);

void
rvec_add(double * a,double * b,double * c);

void
copy_rvec(double * a,double * b);

double norm2(double * a);

void svmul(double a,double *  v1,double *  v2);

void unitv(double *  src,double *  dest);

void oprod(double *  a,double *  b,double *  c);

double iprod(double *  a,double *  b);

double 
cos_angle_no_table(double * a,double * b);

void 
rotvector_from_matrix(double rotmat[3][3], double *v, double *theta);

void 
rotmatrix_from_vector(double * v, double theta, double rotmat[3][3]);

void 
mvmul(double a[3][3],double * src,double * dest);

void 
calc_fit3(double * a1, double * a2, double * a3, 
	 double * b1, double * b2, double * b3,
	 double * cogA, double * cogB, 
	 double * rotvec, double *rotangle);
 


#endif
