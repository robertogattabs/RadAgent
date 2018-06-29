/* Function for calculatinge mean and median dvh */
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <time.h>

/* comparison function for qsort */
int cmpfunc (const void * a, const void * b) {
   const double *da = (const double *) a;
   const double *db = (const double *) b;
   return (*da > *db) - (*da < *db);
}

void meanmediandvh(double *dvh, int *Nrow, int *Nboot, double *meanDvh, 
                   double *sampledvh, int *Nhist, double *medianDvh, double *stepMedian) {
	int n, m, p; /* counters */
  int dvhindex;  /* index of histograms */    
  /* initialize random seed */
  srand ( time(NULL) );  
	for (n=0; n < *Nboot; n++) { /* iteration among bootstrap dvh */
    /* creation of resampled histogram series */    
		for (m=0; m < *Nhist; m++) {  /* iteration among number of histograms */
      dvhindex = rand() % *Nhist; /* create the index of the sampled dvh */
      for (p=0; p < *Nrow; p++) { /* start creating sampled dvh */
        sampledvh[m * *Nrow + p] = dvh[dvhindex * *Nrow + p]; /* set element in the dvh */        
      }      
		}  
    /* calculate mean */
    for (m=0; m < *Nhist; m++) {
      for (p=0; p < *Nrow; p++) {
        meanDvh[n * *Nrow + p] += sampledvh[m * *Nrow + p];        
      }            
    }
    for (p=0; p < *Nrow; p++) {
      meanDvh[n * *Nrow + p] = meanDvh[n * *Nrow + p] / *Nhist;
    }
    /* calculate median */ 
    for (p=0; p < *Nrow; p++) {      // loop along rows
      for (m=0; m < *Nhist; m++) {   // loop along DVHs
        stepMedian[m] = sampledvh[m * *Nrow + p];  // create array of one row of sampled DVH
      }
      /* sort the array */
      qsort(stepMedian, *Nhist, sizeof(double), cmpfunc);
      if(*Nhist % 2 == 0) {
        /* if there is an even number of elements, return mean of the two elements in the middle */
        medianDvh[n * *Nrow + p] = ((stepMedian[*Nhist/2] + stepMedian[*Nhist/2 - 1]) / 2.0);
      } else {
        /* else return the element in the middle */
        medianDvh[n * *Nrow + p] = stepMedian[*Nhist/2];
      }
    }
  }  
}

