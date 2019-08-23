#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "geometry.h"

#define NATYPE 7
static char 
atomtypename [NATYPE] = "HCNOSPX";

static double 
bondlength [NATYPE][NATYPE] = {
  {  0.00, 1.08, 1.05, 1.00, 1.00, 1.00, 1.00 },
  {  1.08, 1.50, 1.40, 1.40, 1.60, 1.60, 1.40 },
  {  1.05, 1.40, 1.40, 1.40, 1.60, 1.60, 1.40 },
  {  1.00, 1.40, 1.40, 1.40, 1.70, 1.70, 1.40 },
  {  1.00, 1.60, 1.60, 1.70, 2.00, 2.00, 1.70 },
  {  1.00, 1.60, 1.60, 1.70, 2.00, 2.00, 1.70 },
  {  1.00, 1.40, 1.40, 1.40, 1.70, 1.70, 1.70 }
};

typedef struct {
  int     nbonds;
  int     bond[10];
  int     branch;
  int     level;
  int     ref0;
  int     ref1;
  int     ref2;
  int     used;
} bondgraph_t;


static void
search_recursively(bondgraph_t *bgraph, 
		   int i, 
		   int ref0, 
		   int ref1, 
		   int ref2, 
		   int level, 
		   int branch)
{
  int j;
  
  if(bgraph[i].branch!=-1 || bgraph[i].used==1)
    return;

  bgraph[i].branch=branch;
  bgraph[i].level=level;
  bgraph[i].ref0=ref0;
  bgraph[i].ref1=ref1;
  bgraph[i].ref2=ref2;

  level++;

  for(j=0;j<bgraph[i].nbonds;j++)
    search_recursively(bgraph,bgraph[i].bond[j],ref1,ref2,i,level,branch);
}


static
void
find_longest_new_path(bondgraph_t *       bgraph, 
		      zmatrix_graph_t *   zgraph,
		      int                 natoms,
		      int *               nfound)
{
  int i,j,k,l;
  int maxlevel,maxidx;
  int nbranch;
  int r0,r1,r2;

  /* reset all unused atoms to branch -1 */
  for(i=0;i<natoms;i++) 
    if(!bgraph[i].used)
      bgraph[i].branch=-1;

  /* start indexing new branches */
  nbranch = zgraph->nbranch;

  if(nbranch==0) {
    /* start from scratch. Never mind reference atoms.
     *
     * Normally all atoms are listed more or less in order in a PDB file,
     * but there is an exception: NH2/CH3 groups at the start. In this
     * case we should start on one of the hydrogens instead.
     */
    if(bgraph[0].bond[0]==1 && bgraph[1].bond[0]==0)
      search_recursively(bgraph,1,-3,-2,-1,1,nbranch);
  
    for(i=0;i<natoms;i++) {
      search_recursively(bgraph,i,-3,-2,-1,1,nbranch);
    }
  } else {
    /* start from previous branches */
    for(i=0;i<natoms;i++) {
      if(bgraph[i].used) {
	for(j=0;j<bgraph[i].nbonds;j++) {
	  k=bgraph[i].bond[j];
	  if(!bgraph[k].used) {
	    /* find references? */
	    r0=bgraph[k].ref0;
	    r1=bgraph[k].ref1;
	    r2=bgraph[k].ref2;

	    /* check if r0/r1 are real atoms (we just started from r2)
	     * otherwise we have to go forward in that branch instead 
	     */
	    if(r0<0) {
	      l=0;
	      while(zgraph->atom[bgraph[i].branch][l] != i)
		l++;
	      r1=zgraph->atom[bgraph[i].branch][l+1];
	      r0=zgraph->atom[bgraph[i].branch][l+2];
	    }
	    search_recursively(bgraph,k,r0,r1,r2,1,nbranch);
	  }
	}
      }
    }
  }

  /* Got several new branches.
   * Find max level and index of it.
   */
  maxlevel=-1;
  for(i=0;i<natoms;i++) 
    if(bgraph[i].branch==nbranch && bgraph[i].level>maxlevel) {
      maxlevel=bgraph[i].level;
      maxidx=i;
    }
  if(maxlevel==-1) {
    fprintf(stderr,"Fatal error. Could only connect %d out of %d atoms.\n",*nfound,natoms);
    exit(1);
  }

  /* Copy the path to the zmatrix graph, and set
   * the branch number in the bond graph.
   */
  zgraph->branchlength=(int *)realloc(zgraph->branchlength,
				      sizeof(int)*(nbranch+1));
  zgraph->atom=(int **)realloc(zgraph->atom,
			       sizeof(int *)*(nbranch+1));
  zgraph->refatom=(i3 **)realloc(zgraph->refatom,
				   sizeof(i3 *)*(nbranch+1));

  zgraph->branchlength[nbranch]=maxlevel;
  
  zgraph->atom[nbranch]=(int *)malloc(sizeof(int)*(maxlevel));
  zgraph->refatom[nbranch]=(i3 *)malloc(sizeof(i3)*(maxlevel));

  /* fill in atoms */
  k=maxidx;
  for(i=maxlevel-1;i>=0;i--) {
    zgraph->atom[nbranch][i]=k;
    zgraph->refatom[nbranch][i][0]=bgraph[k].ref0;
    zgraph->refatom[nbranch][i][1]=bgraph[k].ref1;
    zgraph->refatom[nbranch][i][2]=bgraph[k].ref2;
    bgraph[k].used=1;
    k=bgraph[k].ref2;
  }
  zgraph->nbranch=nbranch+1;

  *nfound += maxlevel;
}


void
create_zmatrix_graph(double *x, pdbinfo_t *pdbinfo, int natoms,
                     zmatrix_graph_t **, int *ngraph)
{
  int i,j;
  int nfound;
  bondgraph_t *bgraph;
  zmatrix_graph_t *zgraph;
  double r,rmax;
  char ch;
  int *atomtype;

  bgraph=(bondgraph_t *)malloc(sizeof(bondgraph_t)*natoms);
  zgraph=(zmatrix_graph_t *)malloc(sizeof(zmatrix_graph_t)*natoms);
  atomtype=(int *)malloc(sizeof(int)*natoms);

  zgraph->atomtobranch=(int *)malloc(sizeof(int)*natoms);

  /* set atom names */
  for(i=0;i<natoms;i++) {
    ch=pdbinfo[i].atomn[0];
    if(ch>='0' && ch<='9')
      ch=pdbinfo[i].atomn[1];
    j=0;
    while(j<NATYPE-1 && ch!=atomtypename[j])
      j++;
    atomtype[i]=j;
    /* if it didn't match any atom it is now equivalent to 'X' */
  }

  for(i=0;i<natoms;i++) 
    bgraph[i].nbonds=0;
  
  /* find all bonds. Yeah, I know this is a horrible N^2 loop... */
  for(i=0;i<natoms;i++) {
    for(j=i+1;j<natoms;j++) {
      r= calc_dist(x,i,j);
      rmax = BOND_MARGIN * bondlength[atomtype[i]][atomtype[j]];
      if(r<rmax) {
	bgraph[i].bond[bgraph[i].nbonds]=j;
	bgraph[i].nbonds++;
	bgraph[j].bond[bgraph[j].nbonds]=i;
	bgraph[j].nbonds++;
      }
    }
  }

  /* Bonds determined. Time to find paths */

  /* Clear all atoms - set 'branch' to -1 */
  for(i=0;i<natoms;i++) 
    bgraph[i].branch=-1;

  for(i=0;i<natoms;i++) 
    bgraph[i].used=0;

  zgraph->nbranch=0;
  zgraph->branchlength=NULL;
  zgraph->atom=NULL;
  zgraph->refatom=NULL;

  nfound=0;
  /* While free atoms remain:
   * 1. Find the longest branch 
   * 2. Extract it, set branch of those atoms to the branch number
   * 3. Repeat.
   */ 
   while(nfound<natoms) 
     find_longest_new_path(bgraph,zgraph,natoms,&nfound);

  for(i=0;i<natoms;i++) 
    zgraph->atomtobranch[i]=bgraph[i].branch;
   
   return zgraph;
}

