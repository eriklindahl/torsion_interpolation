
double
calcrmsd(double *     x1,
     double *     x2,
     int          natoms);

double
calcrmsdindex(double *     x1,
	      double *     x2,
	      int *index,
	      int          nindex);

     
void
fitrmsd(double *  ref_x,
	double *  x,
	int       natoms);

void
fitrmsdindex(double *  ref_x,
	     double *  x,
	     int *     index,
	     int       nindex);
