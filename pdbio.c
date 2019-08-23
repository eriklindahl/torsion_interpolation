#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pdbio.h"

int
pdb_read_first_model(FILE *         fp,
		     double **      x, 
		     pdbinfo_t **   pdbinfo,
		     int *          natoms)
{
    char buffer[255];
    int nat;
    int nalloc;
    
    nat = 0;
    nalloc = 1000;
    *x = malloc(3*sizeof(double)*nalloc);
    *pdbinfo = malloc(sizeof(pdbinfo_t)*nalloc);
    
    while(fgets(buffer,255,fp) && strncmp(buffer,"ENDMDL",6))
    {
        /* Only look at ATOM records, and skip lines with the
         * alternate locator flag in column 17 set */
        if(!strncmp(buffer,"ATOM",4) &&  buffer[16]==' ')
        {
            if(nat==nalloc) {
                nalloc += 1000;
                *x = realloc(*x,3*sizeof(double)*nalloc);
                *pdbinfo = realloc(*pdbinfo,sizeof(pdbinfo_t)*nalloc);
            }
            sscanf(buffer+12,"%5s",((*pdbinfo)[nat].atomn));
            (*pdbinfo)[nat].atomn[4]=0; /* make sure string is terminated */
            
            /* allocate more memory if necessary */
            (*pdbinfo)[nat].chain=buffer[21];
            
            sscanf(buffer+17,"%3s",((*pdbinfo)[nat].resn));
            (*pdbinfo)[nat].resn[4]=0; /* make sure string is terminated */
            sscanf(buffer+22,"%d",&((*pdbinfo)[nat].resi));
            sscanf(buffer+30,"%lf %lf %lf",*x+3*nat,*x+3*nat+1,*x+3*nat+2);
            nat++;
        }
    }
    *natoms = nat;
    return nat;
}



int
pdb_read_next_model(FILE *         fp,
		    double *       x, 
		    pdbinfo_t *    pdbinfo,
		    int            natoms)
{
  char buffer[255];
  int nat;

  nat = 0;

  while(fgets(buffer,255,fp) && nat<natoms) {

    if(!strncmp(buffer,"ENDMDL",6) || !strncmp(buffer,"TER",3)) {
      printf("Fatal error. Your trajectory is broken.\n");
      exit(1);
    }
    
    /* Only look at ATOM records, and skip lines with the 
     * alternate locator flag in column 17 set */
    if(!strncmp(buffer,"ATOM",4) &&  buffer[16]==' ') {
      sscanf(buffer+12,"%5s",(pdbinfo[nat].atomn));
      pdbinfo[nat].atomn[4]=0; /* make sure string is terminated */
       
      pdbinfo[nat].chain=buffer[21];
      
      sscanf(buffer+17,"%3s",(pdbinfo[nat].resn));
      pdbinfo[nat].resn[4]=0; /* make sure string is terminated */
      sscanf(buffer+22,"%d",&(pdbinfo[nat].resi));
      sscanf(buffer+30,"%lf %lf %lf",x+3*nat,x+3*nat+1,x+3*nat+2);
      nat++;
    }
  }
  /* Got the atoms we were looking for. Check that there are no 
   * other atoms before TER/ENDMDL.
   */
  while(fgets(buffer,255,fp) && strncmp(buffer,"ENDMDL",6)) {
    if(!strncmp(buffer,"ATOM",4)) {
      printf("Fatal error. Your trajectory is broken.\n");
      exit(1);
    }
  }
  return nat;
}


int
pdb_write_model(FILE *         fp,
		double *       x, 
		pdbinfo_t *    pdbinfo,
		int            natoms)
{
  int i;
  char lastchain;
  char atomname[5];

  if(natoms>0)
    lastchain=pdbinfo[0].chain;
  else
    lastchain='#'; /* anything but A-Z */

  fprintf(fp,"HEADER Frame generated by pdbio.c\n");
  
  for(i=0;i<natoms;i++) {
    if(pdbinfo[i].chain!=lastchain)
      fprintf(fp,"TER\n");
    if(strlen(pdbinfo[i].atomn)>3)
      sprintf(atomname,"%4s",pdbinfo[i].atomn);
    else
      sprintf(atomname," %-3s",pdbinfo[i].atomn);
      
    fprintf(fp,"ATOM  %5d %4s %3s %c %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
	    i+1,atomname,
	    pdbinfo[i].resn,pdbinfo[i].chain,pdbinfo[i].resi,
	    x[3*i],x[3*i+1],x[3*i+2]);
    lastchain=pdbinfo[i].chain;
  }
  
  fprintf(fp,"TER\nENDMDL\n");

  return natoms;
}

