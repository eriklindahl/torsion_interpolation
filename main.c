#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "pdbio.h"
#include "interpolation.h"
#include "graph.h"
#include "rmsd.h"
#include "spline.h"

int
main(int argc, char **argv)
{
  FILE *fpin,*fpout;
  int i,j,k,n,ret;
  int natoms;
  pdbinfo_t *pdbinfo;
  double **x,*newx,*lastx;
  zdata_t  ****zdata;
    /* zdata[i][j][k][l] refers to:
     * Frame i
     * Molecule j
     * 
     */
  zmatrix_graph_t *graph;
  double ***phi_values;
  double ***phi_der2;
  double lambda;
    int ngraph;
  int frame=0;
  int cnt;
  FILE *fp_phi;
  
  if(argc<4) {
    fprintf(stderr,"Usage:\n\ninterpolate <infile.pdb> <outfile.pdb> <n>\n");
    exit(1);
  }

  fpin=fopen(argv[1],"r");
  fpout=fopen(argv[2],"w");
  sscanf(argv[3],"%d",&n);

  x=(double **)malloc(sizeof(double *)*2);
  fprintf(stderr,"\rReading frame %4d... ",frame++);
  ret=pdb_read_first_model(fpin,&(x[0]),&pdbinfo,&natoms);
    
  if(ret<=0) {
    fprintf(stderr,"Error reading first pdb frame.\n");
    exit(1);
  }

  x[1]=(double *)malloc(sizeof(double)*3*natoms);
  newx=(double *)malloc(sizeof(double)*3*natoms);
  lastx=(double *)malloc(sizeof(double)*3*natoms);

  // create dihedral graph
  create_zmatrix_graph(x[0],pdbinfo,natoms,&graph,&ngraph);
    printf("A3\n");

  /* Allocate a list of zmatrix data, one per molecule */
  zdata=(zdata_t ****)malloc(sizeof(zdata_t ***)*ngraph);
    
  for(i=0;i<ngraph;i++)
  {
      zdata[i][0]=init_zdata(graph[i]);
      calculate_zdata(graph[i],x[0],zdata[i][0]);
  }

    printf("A4\n");

    // read all remaining frames
    frame=1;
    while(pdb_read_next_model(fpin,x[frame],pdbinfo,natoms)>0) {
    fprintf(stderr,"\rReading frame %4d... ",frame);

        
    zdata=(zdata_t ***)realloc(zdata,sizeof(zdata_t **)*(frame+1));
    zdata[frame]=init_zdata(graph);

    // calculate zmatrix data
    calculate_zdata(graph,x[frame],zdata[frame]);

    frame++;

    x=(double **)realloc(x,sizeof(double *)*(frame+1));
    x[frame]=(double *)malloc(sizeof(double)*3*natoms);
  }

  // have the original data.
  fclose(fpin);
  fprintf(stderr,"\nFixing torsions (might take a while)...\n");

  // fix the zmatrix so the interpolation will work
  fix_zdata(graph,x,zdata,natoms,frame);
  printf("fixed\n");
  
  // prepare for interpolation. We need to transpose data first...
  // branch * branchlength * nframes
  phi_values=(double ***)malloc(sizeof(double **)*graph->nbranch);
  phi_der2=(double ***)malloc(sizeof(double **)*graph->nbranch);
  for(i=0;i<graph->nbranch;i++) {
    phi_values[i]=(double **)malloc(sizeof(double *)*graph->branchlength[i]);
    phi_der2[i]=(double **)malloc(sizeof(double *)*graph->branchlength[i]);
    for(j=0;j<graph->branchlength[i];j++) {
      phi_values[i][j]=(double *)malloc(sizeof(double )*(frame));
      phi_der2[i][j]=(double *)malloc(sizeof(double )*(frame));
    }
  }
  // move data and create splines
  for(i=0;i<graph->nbranch;i++) {
      for(j=0;j<graph->branchlength[i];j++) {
          for(k=0;k<frame;k++)
              phi_values[i][j][k]=zdata[k][i][j].phi;
      create_spline(phi_values[i][j],frame,1e30,1e30,phi_der2[i][j]);
    }
  }
  
  // write first frame 
  pdb_write_model(fpout,x[0],pdbinfo,natoms);
  memcpy(lastx,x[0],sizeof(double)*3*natoms);
  fprintf(stderr,"\rWriting frame 0... ");
  
  for(i=1;i<frame;i++) {
      
    for(j=1;j<n;j++) {

      lambda=((double)j)/((double)n);
      interpolate_with_splines(graph,zdata[i-1],zdata[i],phi_values,phi_der2,lambda,i-1.0,frame,x[i-1],x[i],newx);
      fitrmsd(lastx,newx,natoms);
      fprintf(stderr,"\rWriting original frame %4d, interpolation %4d / %4d... ",i-1,j,n);
      pdb_write_model(fpout,newx,pdbinfo,natoms);
      memcpy(lastx,newx,sizeof(double)*3*natoms);
    }
    fitrmsd(lastx,x[i],natoms);
    fprintf(stderr,"\rWriting original frame %4d, interpolation %4d / %4d... ",i-1,n,n);
    pdb_write_model(fpout,x[i],pdbinfo,natoms);
    
    memcpy(lastx,x[i],sizeof(double)*3*natoms);
  }
  fclose(fpout);
  
  fprintf(stderr,"\nDone.\n");

  return 0;
}

