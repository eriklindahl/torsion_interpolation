
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "graph.h"
#include "interpolation.h"
#include "geometry.h"
#include "rmsd.h"
#include "spline.h"

zdata_t **
init_zdata(zmatrix_graph_t *graph)
{
  int i;
  zdata_t **p;

  p=(zdata_t **)malloc(sizeof(zdata_t *)*graph->nbranch);
  
  for(i=0;i<graph->nbranch;i++)
    p[i]=(zdata_t *)malloc(sizeof(zdata_t)*graph->branchlength[i]);

  return p;
}



/* Before calling this routine, 'zdata' must be allocated
 * as a pointer to a list of length nbranch, where each entry
 * is a pointer to a list of length branchlength[i].
 * Use init_zdata() for simplicity.
 */
void
calculate_zdata(zmatrix_graph_t *   graph,
		double *            x,
		zdata_t **          zdata)
{
  int i,j,at,r0,r1,r2;
  double r,theta,phi;

  for(i=0;i<graph->nbranch;i++) 
    for(j=0;j<graph->branchlength[i];j++) {
      at = graph->atom[i][j];
      r0 = graph->refatom[i][j][0];
      r1 = graph->refatom[i][j][1];
      r2 = graph->refatom[i][j][2];
      
      if(r0<0) {
	r=theta=phi=0;
      } else {
	r = calc_dist(x,at,r2);
	theta = calc_angle(x,at,r2,r1);
	phi = calc_torsion(x,at,r2,r1,r0);
      }
      zdata[i][j].r = r;
      zdata[i][j].theta = theta;
      zdata[i][j].phi = phi;
    } 
}


static void
zmatrix_to_cartesian(int         length,
		     int *       atom,
		     i3 *        refatom,
		     zdata_t *   zdataA,
		     zdata_t *   zdataB,
		     double      lambda,
		     double *    x)
{
  int a,r2,r1,r0;
  double u[3],v[3],tmpx[3],ux[3],w[3];
  int i;
  double r,theta,phi,scratch;

  for(i=0;i<length;i++) {
    a=atom[i];
    r0=refatom[i][0];
    r1=refatom[i][1];
    r2=refatom[i][2];
    
    if(r0<0)
      continue;
    
    r=(1.0-lambda)*zdataA[i].r + lambda*zdataB[i].r;
    theta=(1.0-lambda)*zdataA[i].theta + lambda*zdataB[i].theta;
    phi= (1.0-lambda)*zdataA[i].phi + lambda*zdataB[i].phi;

    rvec_sub(x+3*r2,x+3*r1,u);
    rvec_sub(x+3*r0,x+3*r1,v);
    
    unitv(u,u);
    scratch=iprod(u,v);
    svmul(scratch,u,ux);
    rvec_sub(v,ux,v);
    unitv(v,v);
    oprod(u,v,w);
    scratch=r*sin(theta);
    
    svmul(-r*cos(theta),u,u);
    svmul(scratch*cos(phi),v,v);
    svmul(scratch*sin(phi),w,w);
    
    rvec_add(x+3*r2,u,tmpx);
    rvec_add(tmpx,v,tmpx);
    rvec_add(tmpx,w,tmpx);
    
    copy_rvec(tmpx,x+3*a);
  }
}


static void
zmatrix_to_cartesian_with_splines(int         length,
				  int *       atom,
				  i3 *        refatom,
				  zdata_t *   zdataA,
				  zdata_t *   zdataB,
				  double **   phi_values,
				  double **   phi_der2,
				  double      lambda,
				  int         wholestep,
				  int         nframe,
				  double *    x)
{
  int a,r2,r1,r0;
  double u[3],v[3],tmpx[3],ux[3],w[3];
  int i;
  double splinex;
  double r,theta,phi,scratch;

  for(i=0;i<length;i++) {
    a=atom[i];
    r0=refatom[i][0];
    r1=refatom[i][1];
    r2=refatom[i][2];
    
    if(r0<0)
      continue;
    
    r=(1.0-lambda)*zdataA[i].r + lambda*zdataB[i].r;
    theta=(1.0-lambda)*zdataA[i].theta + lambda*zdataB[i].theta;

    /* get phi from spline interpolation */
    splinex=(double)wholestep+lambda;
    phi=spline_interpolation(phi_values[i],phi_der2[i],nframe,splinex);

    rvec_sub(x+3*r2,x+3*r1,u);
    rvec_sub(x+3*r0,x+3*r1,v);
    
    unitv(u,u);
    scratch=iprod(u,v);
    svmul(scratch,u,ux);
    rvec_sub(v,ux,v);
    unitv(v,v);
    oprod(u,v,w);
    scratch=r*sin(theta);
    
    svmul(-r*cos(theta),u,u);
    svmul(scratch*cos(phi),v,v);
    svmul(scratch*sin(phi),w,w);
    
    rvec_add(x+3*r2,u,tmpx);
    rvec_add(tmpx,v,tmpx);
    rvec_add(tmpx,w,tmpx);
    
    copy_rvec(tmpx,x+3*a);
  }
}



void 
interpolate_linear  (zmatrix_graph_t *   graph,
		     zdata_t **          zdataA,
		     zdata_t **          zdataB,
		     double              lambda,
		     double *            xA,
		     double *            xB,
		     double *            x)
{
  int i;
  int r0,r1,r2;
  double rotangle,rotvec[3];
  double cogA[3],cogB[3];
  double tv1[3],tv2[3],tv3[3],tv4[3],tv5[3],tv6[3],newcom[3];
  double rotmat[3][3];
  
  /* First interpolate the first three reference atoms 
   * in the first branch using rigid rotation.
   */

  r0=graph->atom[0][0];
  r1=graph->atom[0][1];
  r2=graph->atom[0][2];

  calc_fit3(xA+3*r0,xA+3*r1,xA+3*r2,
	    xB+3*r0,xB+3*r1,xB+3*r2,
	    cogA,cogB,rotvec,&rotangle);

  svmul(1.0-lambda,cogA,tv1);
  svmul(lambda,cogB,tv2);
  rvec_add(tv1,tv2,newcom);
  rotmatrix_from_vector(rotvec,lambda*rotangle,rotmat);
  rvec_sub(xA+3*r0,cogA,tv1);
  rvec_sub(xA+3*r1,cogA,tv2);
  rvec_sub(xA+3*r2,cogA,tv3);    
  mvmul(rotmat,tv1,tv4);
  mvmul(rotmat,tv2,tv5);
  mvmul(rotmat,tv3,tv6);
  rvec_add(newcom,tv4,x+3*r0);
  rvec_add(newcom,tv5,x+3*r1);
  rvec_add(newcom,tv6,x+3*r2);
 
 
  /* Now do torsion interpolation in the branches.
   */

  for(i=0;i<graph->nbranch;i++) {
    zmatrix_to_cartesian(graph->branchlength[i],
			 graph->atom[i],graph->refatom[i],zdataA[i],zdataB[i],lambda,x);
  }
}




void 
interpolate_with_splines  (zmatrix_graph_t *   graph,
			   zdata_t **          zdataA,
			   zdata_t **          zdataB,
			   double ***          phi_values,
			   double ***          phi_der2,
			   double              lambda,
			   int                 wholeframe,
			   int                 nframe,
			   double *            xA,
			   double *            xB,
			   double *            x)
{
  int i;
  int r0,r1,r2;
  double rotangle,rotvec[3];
  double cogA[3],cogB[3];
  double tv1[3],tv2[3],tv3[3],tv4[3],tv5[3],tv6[3],newcom[3];
  double rotmat[3][3];
  
  /* First interpolate the first three reference atoms 
   * in the first branch using rigid rotation.
   */

  r0=graph->atom[0][0];
  r1=graph->atom[0][1];
  r2=graph->atom[0][2];

  calc_fit3(xA+3*r0,xA+3*r1,xA+3*r2,
	    xB+3*r0,xB+3*r1,xB+3*r2,
	    cogA,cogB,rotvec,&rotangle);

  svmul(1.0-lambda,cogA,tv1);
  svmul(lambda,cogB,tv2);
  rvec_add(tv1,tv2,newcom);
  rotmatrix_from_vector(rotvec,lambda*rotangle,rotmat);
  rvec_sub(xA+3*r0,cogA,tv1);
  rvec_sub(xA+3*r1,cogA,tv2);
  rvec_sub(xA+3*r2,cogA,tv3);    
  mvmul(rotmat,tv1,tv4);
  mvmul(rotmat,tv2,tv5);
  mvmul(rotmat,tv3,tv6);
  rvec_add(newcom,tv4,x+3*r0);
  rvec_add(newcom,tv5,x+3*r1);
  rvec_add(newcom,tv6,x+3*r2);
 
 
  /* Now do torsion interpolation in the branches.
   */

  for(i=0;i<graph->nbranch;i++) {
    zmatrix_to_cartesian_with_splines(graph->branchlength[i],
				      graph->atom[i],graph->refatom[i],zdataA[i],zdataB[i],
				      phi_values[i],phi_der2[i],lambda,wholeframe,nframe,x);
  }
}



void 
copy_zdata       (zmatrix_graph_t *    graph,
		  zdata_t **           src,
		  zdata_t **           dst)
{
  int i,j;

  for(i=0;i<graph->nbranch;i++) {
    for(j=0;j<graph->branchlength[i];j++) {
      dst[i][j].r=src[i][j].r;
      dst[i][j].theta=src[i][j].theta;
      dst[i][j].phi=src[i][j].phi;
    }
  }
}

    
void
fix_zdata(zmatrix_graph_t *graph, double **x, zdata_t ***zdata, int natoms, int nframe)
{
  int i,j,k,n,l,frame;
  int pos;
  double *newx;
  double *xinter[11];
  double dphi;
  int ntest,ncomb,comb;
  int * testpos=NULL;
  double * testval=NULL;
  double * bestval=NULL;
  double rmsd,lowest_sum_crms;
  double acc_dphi;
  int startpos;
  double crms,sum_crms;
  double otherdphi;
  int branchatom,otherbranch;
  int *tmpidx;
  int nidx;
  
  for(i=0;i<11;i++) 
    xinter[i]=(double *)malloc(sizeof(double)*3*natoms);

  newx=(double *)malloc(sizeof(double)*3*natoms);
  tmpidx=(int *)malloc(sizeof(int)*natoms);
  
  testpos=(int *)malloc(sizeof(int)*10000);
  testval=(double *)malloc(sizeof(int)*10000);
  bestval=(double *)malloc(sizeof(int)*10000);

  for(frame=1;frame<nframe;frame++) {
    
      /* Check each z-matrix branch separately */
      for(i=0;i<graph->nbranch;i++) {
          
          ntest=0;
          ncomb=1;
          acc_dphi=0;
          
          /* New branch. */
          if(i>0) {
              /* Preserve stereochemistry around branch point if this
               * is a branch starting from an earlier branch.
               */
              branchatom=graph->refatom[i][0][2];
              /* find our branch sibling - i.e. the other branch from this point */
              otherbranch=graph->atomtobranch[branchatom];
              j=0;
              while(j<graph->branchlength[otherbranch] &&
                    graph->refatom[otherbranch][j][2]!=branchatom)
                  j++;
              
              if(graph->refatom[otherbranch][j][0]<0) {
                  /* beginning of first chain. We can find
                  * a useable torsion to compare with by forwarding one step.
                  */
                  j++;
              }
              otherdphi=zdata[frame][otherbranch][j].phi-zdata[frame-1][otherbranch][j].phi;
              
              dphi=zdata[frame][i][0].phi-zdata[frame-1][i][0].phi;
           
              /* Remove basic periodicity of torsion angle */
              while( dphi-otherdphi > M_PI) 
                  dphi -= 2*M_PI;
              while( dphi-otherdphi < -M_PI)
                  dphi += 2*M_PI;
              
              zdata[frame][i][0].phi=zdata[frame-1][i][0].phi+dphi;
          }
          
          if(i==0)
              startpos=0;
          else
              startpos=1;
          
          for(j=startpos;j<graph->branchlength[i];j++) {
              dphi=zdata[frame][i][j].phi-zdata[frame-1][i][j].phi;
              
              /* Reset basic periodicity of torsions */
              while(dphi>M_PI)
                  dphi -= 2*M_PI;
              while(dphi<-M_PI)
                  dphi += 2*M_PI;

              zdata[frame][i][j].phi=zdata[frame-1][i][j].phi+dphi;
              /* We mark all rotations that are 'big', to make sure
               * they are not interpolated in the wrong direction.
               * For now we consider everything larger than PI/2 as big...
               */
              if(fabs(dphi)>CHECK_TORSION_LIMIT) {
                  testpos[ntest]=j;
                  testval[ntest]=dphi;
                  ntest++;
                  ncomb = ncomb * 2;
              }
          }
          /* If there is more than one big rotation we need to check combinations.
           * This is done by trying all alternatives, and calculating the sum of 
           * cRMS between interpolated structures between the two original states.
           * The lowest sum corresponds to the best interpolation directions.    
           */
          if(ntest>0) {
              /* Find best direction (CW/CCW) for large torsion motions */
              lowest_sum_crms=1e10; 
              for(comb=0;comb<ncomb;comb++) {
                  /* set values */
                  for(n=0;n<ntest;n++) 
                  {
                      pos    = testpos[n];
                      dphi   = testval[n];
                      
                      if(comb & (1 << n)) 
                      {
                          /* flag is 1, this dphi should be positive */
                          if(dphi<0)
                              dphi += 2*M_PI;
                      }
                      else
                      {
                          /* this dphi should be negative */
                          if(dphi>0)
                              dphi -= 2*M_PI;
                      }
                      testval[n]=dphi;
                      zdata[frame][i][pos].phi=zdata[frame-1][i][pos].phi+dphi;
                  }
                  
	  
                  // interpolate 9 new positions between states
                  for(l=1;l<10;l++)
                      interpolate_linear(graph,zdata[frame-1],zdata[frame],0.1*l,
                                         x[frame-1],x[frame],xinter[l]);

                  // create a temporary index
                  for(l=0;l<graph->branchlength[i];l++)
                      tmpidx[l] = graph->atom[i][l];
                  
                  nidx = graph->branchlength[i];
                  
                  if(i>0) {
                      /* Add our reference atoms (nonexistant for 1st branch) */
                      tmpidx[nidx++] = graph->refatom[i][0][0];
                      tmpidx[nidx++] = graph->refatom[i][0][1];
                      tmpidx[nidx++] = graph->refatom[i][0][2];
                   }

                  /* Copy original coordinates to xinter array */
                  memcpy(xinter[0],x[frame-1],sizeof(double)*3*natoms);
                  memcpy(xinter[10],x[frame],sizeof(double)*3*natoms);
                  
                         
                  fitrmsdindex(xinter[0],xinter[10],tmpidx,nidx);
                  crms = calcrmsdindex(xinter[0],xinter[10],tmpidx,nidx);
                  
                  sum_crms = 0;
                  
                  for(l=0;l<10;l++)
                  {
                      fitrmsdindex(xinter[l],xinter[l+1],tmpidx,nidx);
                      sum_crms += calcrmsdindex(xinter[l],xinter[l+1],tmpidx,nidx);
                  } 
                  
                  // sum_crms = calcrmsd(x[frame-1],xinter[1],natoms);
                  
                  // for(l=1;l<9;l++)
                  //    sum_crms += calcrmsd(xinter[l],xinter[l+1],natoms);
                  
                  // sum_crms += calcrmsd(xinter[9],x[frame],natoms);
                  
                  if(sum_crms<lowest_sum_crms) {
                      lowest_sum_crms=sum_crms;
                      for(k=0;k<ntest;k++)
                          bestval[k]=testval[k];
                  }
              }
              /* Check if there could be errors.
               *
               * Warn if the cRMS sum is more than three times larger
               * compared to the cRMS between the original frames.
               * Skip warnings for small moves (i.e. sum less than 5AA)    
               */   
              if((lowest_sum_crms>3.0*crms) && lowest_sum_crms>5.0)
              {
                  printf("Warning: Large sum of interpolated cRMS (%.3f) compared to cRMS between\n"
                         "original frames %d and %d (%.3f). This *could* be an error, but usually isn't.\n"
                         "You can check the interpolated trajectory in pymol - our framecount starts on 0.\n",
                         lowest_sum_crms,frame-1,frame,crms);
              }
              // set values to the best values
              for(n=0;n<ntest;n++) {
                  pos    = testpos[n];
                  dphi   = bestval[n];
                  zdata[frame][i][pos].phi=zdata[frame-1][i][pos].phi+dphi;
              }
          }
      }
  }  
  free(testval);
  free(bestval);
  free(testpos);
  free(newx);
  free(tmpidx);
  
  for(i=0;i<10;i++)
      free(xinter[i]);
}
