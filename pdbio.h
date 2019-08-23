#include <stdio.h>

#ifndef _PDBIO_H_
#define _PDBIO_H_

typedef struct {
  int       resi;
  char      resn[5];
  char      atomn[5];
  char      chain;
} pdbinfo_t;


/* read the first frame - i.e. until we find a ENDMDL record */
int
pdb_read_first_model (FILE *         fp,
		      double **      x, 
		      pdbinfo_t **   pdbinfo,
		      int *          natoms);

/* read the next natoms ATOM records, and check that they are
 * followed by ENDMDL
 */
int
pdb_read_next_model  (FILE *         fp,
		      double *       x, 
		      pdbinfo_t *    pdbinfo,
		      int            natoms);

int
pdb_write_model      (FILE *         fp,
		      double *       x, 
		      pdbinfo_t *    pdbinfo,
		      int            natoms);

#endif
