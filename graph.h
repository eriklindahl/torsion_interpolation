#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "pdbio.h"

typedef int i3[3];



/* Allow bond lengths to be 10% too long */
#define BOND_MARGIN 1.1


typedef struct {
  int                  nbranch;    
  int *                branchlength;
  int **               atom;
  i3 **                refatom; /* 3 atoms for each entry in each branch */
  int *                atomtobranch;
} zmatrix_graph_t;



void
create_zmatrix_graph(double *x, pdbinfo_t *pdbinfo,
                     int natoms,
                     zmatrix_graph_t **,
                     int *ngraph);


#endif
