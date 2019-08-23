#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_


/* Torsion angles that change by more than this value
 * are checked to avoid interpolating in the wrong
 * direction. Note that the computational cost goes as
 * 2^(N), where N is the number of torsion angles to check.
 * In other words - if you make this too low it is going to
 * take forever to run the program, and if it is too high
 * you might do bad interpolations. PI/2 is a reasonable
 * initial value, but you can change it.
 */
#define CHECK_TORSION_LIMIT (0.5*M_PI)


#include "graph.h"

typedef struct {
  double               r;
  double               theta;
  double               phi;
} zdata_t;


zdata_t **
init_zdata(zmatrix_graph_t *graph);


/* Before calling this routine, the zdata must be allocated
 * as a pointer to a list of length nbranches, where each entry
 * is a pointer to a list of length branchlength[i].
 * Use init_torsiondata() for simplicity.
 */
void
calculate_zdata(zmatrix_graph_t *graph, double *x, zdata_t **newz);

void
fix_zdata(zmatrix_graph_t *graph, double **x, zdata_t ***zdata, int natoms, int nframe);

void 
interpolate_linear  (zmatrix_graph_t *   graph,
		     zdata_t **          zdataA,
		     zdata_t **          zdataB,
		     double              lambda,
		     double *            xA,
		     double *            xB,
		     double *            x);

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
			   double *            x);

void 
copy_zdata          (zmatrix_graph_t *   graph,
		     zdata_t **          src,
		     zdata_t **          dst);

    

#endif
