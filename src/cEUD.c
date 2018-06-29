/* EUD algorithm written in C */
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>


void cEUD(double *dosebin, double *volumebin, double *afactor, int *Nrow, int *Ncol, double *ceud) {
	int n, m, offset;
	double *punt; 		// pointer to temporary dose bin structure
	double addendum; 	// addendum of dose*volume
	punt=(double *)calloc(*Nrow, sizeof(double)); // allocation of pointer
	for (n=0; n<= *Nrow - 1; n++) { // iteration for pre-calculate dose bins
		punt[n]=pow(dosebin[n], *afactor);		
	}
	for (m=0; m<= *Ncol - 1; m++) {  // iteration for EUD calculation
		offset = *Nrow * m;
		addendum = 0;
		for (n = 0; n<= *Nrow - 1; n++) {
			addendum += punt[n] * (volumebin[n + offset]);
		}
		ceud[m] = pow(addendum, 1/(*afactor));
	}	
	free(punt);
}
