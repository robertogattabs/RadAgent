#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "R.h"
#include <Rinternals.h>
#include <complex.h>
#include <limits.h>
#define min(a, b) (((a) < (b)) ? (a) : (b));

struct pointInSpace{
  double x;
  double y;
  double z;
};
struct DICOM_OrientationMatrix{
  double a11;		// 1st element of ImageOrientationPatient
  double a21;		// 2nd element of ImageOrientationPatient
  double a31;		// 3rd element of ImageOrientationPatient	
  double a12;		// 4th element of ImageOrientationPatient
  double a22;		// 5th element of ImageOrientationPatient
  double a32;		// 6th element of ImageOrientationPatient
  double Sx;		// 1st element of ImagePositionPatient
  double Sy;		// 2nd element of ImagePositionPatient
  double Sz;		// 3nd element of ImagePositionPatient
  double XpixelSpacing;	// 1st element of Pixel Spacing
  double YpixelSpacing;	// 2nd element of Pixel Spacing
  struct DICOM_OrientationMatrix *next;
  struct DICOM_OrientationMatrix *prev;
};
struct _c_data{
    // data passed as parameters
    int xNVoxel,yNVoxel,zNVoxel;
    double xDim,yDim,zDim;
    int newXNVoxel,newYNVoxel,newZNVoxel;
    double newXDim,newYDim,newZDim;
    // Voxel number for each vertex of the recognized Cube
    int xInf,xSup,yInf,ySup,zInf,zSup;
};

/*
void nOrderMoment(double *matrix, int *Nx, int *Ny, int *Nz,double *order,double *mX,double *mY,double *mZ ) {

	mX=0;	mY=0;	mZ=0;
	int x,y,z;
	int numberOfPoints=0;

	for(z=0; z<(*Nz); z++) {
		for(y=0; z<(*Ny); y++) {
			for(x=0; x<(*Nz); x++) {		
				if ( matrix[x + (y * (*Nx)) +  z * ( x + y )] !=0)   {
					*mX += pow(x, *order );
					*mY += pow(y, *order );
					*mZ += pow(z, *order );
					numberOfPoints++;
				}
			}
		}
	}
	*mX = *mX/numberOfPoints;
	*mY = *mY/numberOfPoints;
	*mZ = *mZ/numberOfPoints;
}
 */
//
// _c_TrilinearInterpolation
//
// Effettua l'interpolazione considerando i valori ai vertici e le dimensioni fisiche del cubo
//
double _c_TrilinearInterpolation(double x0y0z0, double x0y0z1, double x0y1z0, double x0y1z1, double x1y0z0, double x1y0z1, double x1y1z0, double x1y1z1, double x0,double y0, double z0, double dx1x0, double dy1y0, double dz1z0, double x, double y,double z) {
  double xd,yd,zd,c00,c01,c10,c11,c0,c1,c;
  xd = (x-x0)/dx1x0;
  yd = (y-y0)/dy1y0;
  zd = (z-z0)/dz1z0;
  c00 = x0y0z0*(1-xd)+x1y0z0*xd;
  c10 = x0y1z0*(1-xd)+x1y1z0*xd;
  c01 = x0y0z1*(1-xd)+x1y0z1*xd;
  c11 = x0y1z1*(1-xd)+x1y1z1*xd;
  
  c0 = c00*(1-yd)+c10*yd;
  c1 = c01*(1-yd)+c11*yd;
  c = c0*(1-zd)+c1*zd;
  return c;
}

//
// _c_getCubeVertex
//
// Prende le coordinate in floating point del punto da mappare e restituisce l'intero
// corrispondente ai centroidi relativi alle coordinate dei vertici del cubo da interpolare
// La formula per calcolare il Kinf. è   Kinf = int( (2*xPos - dx) / (2 * dx) )
//
void _c_getCubeVertex(struct _c_data * punt, double xPos, double yPos, double zPos, int ct) {

    punt->xInf = (int)((2 * xPos - punt->xDim) / ( 2 * punt->xDim));
    punt->yInf = (int)((2 * yPos - punt->yDim) / ( 2 * punt->yDim));
    punt->zInf = (int)((2 * zPos - punt->zDim) / ( 2 * punt->zDim));
    punt->xSup = punt->xInf + 1;
    punt->ySup = punt->yInf + 1;
    punt->zSup = punt->zInf + 1;

    // se xPos < punt->xDim significa che sono nel bordo. Considera allora il valore del primo valore di X per x0 e x1
    // (in sostanza proietta all'esterno il valore dei pixel). Analogo se maggiore.
    // STO ANDANDO IN OVERRIDE RISPETTO A QUANDO EVENTUALMENTE CALCOLATO PRIMA
    if( (2*xPos-punt->xDim) < punt->xDim/2 ) {punt->xInf = 0; punt->xSup = 0;}
    if( (2*yPos-punt->yDim) < punt->yDim/2 ) {punt->yInf = 0; punt->ySup = 0;}
    if( (2*zPos-punt->zDim) < punt->zDim/2  ) {punt->zInf = 0; punt->zSup = 0;}

    // Faccio il reciproco di xPos per vedere se non sto' toccando l'altro estremo
    xPos = punt->xDim * punt->xNVoxel - xPos;
    yPos = punt->yDim * punt->yNVoxel - yPos;
    zPos = punt->zDim * punt->zNVoxel - zPos;
    if( (2*xPos-punt->xDim) < punt->xDim/2 ) {punt->xInf = punt->xNVoxel - 1; punt->xSup = punt->xNVoxel - 1;}
    if( (2*yPos-punt->yDim) < punt->yDim/2 ) {punt->yInf = punt->yNVoxel - 1; punt->ySup = punt->yNVoxel - 1;}
    if( (2*zPos-punt->zDim) < punt->zDim/2  ) {punt->zInf = punt->zNVoxel - 1; punt->zSup = punt->zNVoxel - 1;}
}



void newnewtrilinearInterpolator(
  int *NXold, int *NYold,int *NZold,
  int *NXnew, int *NYnew,int *NZnew,
  double *oldXps, double *oldYps, double *oldZps,
  double *newXps, double *newYps, double *newZps,
  double *values,double *returnMatrix) {
  int ct;
  int zVoxelProgressivo,yVoxelProgressivo,xVoxelProgressivo;
  double zPos,yPos,xPos,valoreCalcolato;
  struct _c_data *punt;
  
  // Alloca un puntatore alla struttura _c_data
  punt = (struct _c_data *)calloc(1,sizeof(struct _c_data));
  if( punt == NULL ) return;
  
  // calcola il nuovo passo e copia i valori in _c_data: l'idea è quella di ridurre il passaggio parametri
  // fra funzioni ed il clone delle variabili per preservare memoria 
  punt->newXDim = *newXps;  punt->newYDim = *newYps;  punt->newZDim = *newZps;
  punt->xNVoxel = *NXold; punt->yNVoxel = *NYold; punt->zNVoxel = *NZold;
  punt->xDim = *oldXps; punt->yDim = *oldYps; punt->zDim = *oldZps;
  int maxNewXVoxel,maxNewYVoxel,maxNewZVoxel;
  
  // Pulisci la matrice di destinazione
  for(int ct=0;ct< ((*NXnew)*(*NYnew)*(*NZnew));ct++ ) returnMatrix[ct]=0;  
  ct=0;
  for( zPos = punt->newZDim/2, zVoxelProgressivo=0; zVoxelProgressivo < *NZnew; zPos+=punt->newZDim, zVoxelProgressivo++ ) {
    for( yPos = punt->newYDim/2, yVoxelProgressivo=0; yVoxelProgressivo < *NYnew; yPos+=punt->newYDim, yVoxelProgressivo++ ) {
      for( xPos = punt->newXDim/2, xVoxelProgressivo=0; xVoxelProgressivo < *NXnew; xPos+=punt->newXDim, xVoxelProgressivo++ ) {
        
        // acquisisci i dati relativi ai voxel della vecchia matrice i cui centroidi sono
        // ai vertici del cubo da interpolare
        _c_getCubeVertex( punt , xPos , yPos , zPos ,ct );
        
        // reperisci i valori di tali voxel dalla vecchia matrice
        // ed effettua l'interpolazione rispetto a tali centroidi
        valoreCalcolato = _c_TrilinearInterpolation(
            values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z0 (sample value)
            values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z1 (sample value)
            values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z0 (sample value)
            values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z1 (sample value)
            values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z0 (sample value)
            values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z1 (sample value)
            values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z0 (sample value)
            values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z1 (sample value)
            punt->xInf*(*oldXps)+(*oldXps)/2, //x0
            punt->yInf*(*oldYps)+(*oldYps)/2, //y0,
            punt->zInf*(*oldZps)+(*oldZps)/2, //z0,
            *oldXps, //dx1x0,
            *oldYps, //dy1y0,
            *oldZps, //dz1z0,
            xPos, yPos, zPos);

        // memorizza il risultato nella nuova struttura
        returnMatrix[ xVoxelProgressivo + yVoxelProgressivo * (*NXnew) + zVoxelProgressivo * ((*NYnew) * (*NXnew))] = valoreCalcolato;        
        // if( valoreCalcolato!=0 ) printf("\nx = %d, y= %d, z= %d",xVoxelProgressivo,yVoxelProgressivo,zVoxelProgressivo);
      }
      if(maxNewYVoxel<yVoxelProgressivo) maxNewYVoxel = yVoxelProgressivo;
  }
    if(maxNewZVoxel<zVoxelProgressivo) maxNewZVoxel = zVoxelProgressivo;
}
  free(punt);
  
  }
/*
void newnewtrilinearInterpolator_onGivenMatrix(
  int *NXold, int *NYold,int *NZold,
  int *NXnew, int *NYnew,int *NZnew,
  double *oldXps, double *oldYps, double *oldZps,
  double *newXps, double *newYps, double *newZps,
  double *values,double *returnMatrix) {
  int ct;
  int zVoxelProgressivo,yVoxelProgressivo,xVoxelProgressivo;
  double zPos,yPos,xPos,valoreCalcolato;
  struct _c_data *punt;
  
  // Alloca un puntatore alla struttura _c_data
  punt = (struct _c_data *)calloc(1,sizeof(struct _c_data));
  if( punt == NULL ) return;
  
  // calcola il nuovo passo e copia i valori in _c_data: l'idea è quella di ridurre il passaggio parametri
  // fra funzioni ed il clone delle variabili per preservare memoria 
  punt->newXDim = *newXps;  punt->newYDim = *newYps;  punt->newZDim = *newZps;
  punt->xNVoxel = *NXold; punt->yNVoxel = *NYold; punt->zNVoxel = *NZold;
  punt->xDim = *oldXps; punt->yDim = *oldYps; punt->zDim = *oldZps;
  int maxNewXVoxel,maxNewYVoxel,maxNewZVoxel;
  
  // Pulisci la matrice di destinazione
  for(int ct=0;ct< ((*NXnew)*(*NYnew)*(*NZnew));ct++ ) returnMatrix[ct]=0;  
  ct=0;
  for( zPos = punt->newZDim/2, zVoxelProgressivo=0; zVoxelProgressivo < *NZnew; zPos+=punt->newZDim, zVoxelProgressivo++ ) {
    for( yPos = punt->newYDim/2, yVoxelProgressivo=0; yVoxelProgressivo < *NYnew; yPos+=punt->newYDim, yVoxelProgressivo++ ) {
      for( xPos = punt->newXDim/2, xVoxelProgressivo=0; xVoxelProgressivo < *NXnew; xPos+=punt->newXDim, xVoxelProgressivo++ ) {
        
        // acquisisci i dati relativi ai voxel della vecchia matrice i cui centroidi sono
        // ai vertici del cubo da interpolare
        _c_getCubeVertex( punt , xPos , yPos , zPos ,ct );
        
        // reperisci i valori di tali voxel dalla vecchia matrice
        // ed effettua l'interpolazione rispetto a tali centroidi
        valoreCalcolato = _c_TrilinearInterpolation(
          values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z1 (sample value)
          punt->xInf*(*oldXps)+(*oldXps)/2, //x0
          punt->yInf*(*oldYps)+(*oldYps)/2, //y0,
          punt->zInf*(*oldZps)+(*oldZps)/2, //z0,
          *oldXps, //dx1x0,
          *oldYps, //dy1y0,
          *oldZps, //dz1z0,
          xPos, yPos, zPos);
        
        // memorizza il risultato nella nuova struttura
        returnMatrix[ xVoxelProgressivo + yVoxelProgressivo * (*NXnew) + zVoxelProgressivo * ((*NYnew) * (*NXnew))] = valoreCalcolato;        
        // if( valoreCalcolato!=0 ) printf("\nx = %d, y= %d, z= %d",xVoxelProgressivo,yVoxelProgressivo,zVoxelProgressivo);
      }
      if(maxNewYVoxel<yVoxelProgressivo) maxNewYVoxel = yVoxelProgressivo;
    }
    if(maxNewZVoxel<zVoxelProgressivo) maxNewZVoxel = zVoxelProgressivo;
  }
  free(punt);
  
}
*/
/*
 * NAME: newnewtrilinearInterpolator_onGivenPoints
 * NXold,NYold,NZold : numero di voxel della matrice in ingresso
 * oldXps,oldYps,oldZps : pixelSpacing lungo le 3 direzioni della matrice in ingresso
 * newXpts,newYpts,newZpts : lista dei punti su cui andare ad interpolare
 * pts_x,pts_y,pts_z : numero di voxel della nuova matrice
 * values : matrice in ingresso (in forma di array)
 * returnMatrix : matrice di uscita (in forma di array)
 */
/*
void newnewtrilinearInterpolator_onGivenPoints(
  int *NXold, int *NYold,int *NZold,
  double *oldXps, double *oldYps, double *oldZps,
  double *newXpts, double *newYpts, double *newZpts, int *pts_x, int *pts_y, int *pts_z,
  double *values,double *returnMatrix) {
  int ct;
  int zVoxelProgressivo,yVoxelProgressivo,xVoxelProgressivo;
  double zPos,yPos,xPos,valoreCalcolato;
  struct _c_data *punt;
  int xIndex,yIndex,zIndex;

  // Alloca un puntatore alla struttura _c_data
  punt = (struct _c_data *)calloc(1,sizeof(struct _c_data));
  if( punt == NULL ) return;
  
  // calcola il nuovo passo e copia i valori in _c_data: l'idea è quella di ridurre il passaggio parametri
  // fra funzioni ed il clone delle variabili per preservare memoria 
  //punt->newXDim = *newXps;  punt->newYDim = *newYps;  punt->newZDim = *newZps;
  punt->xNVoxel = *NXold; punt->yNVoxel = *NYold; punt->zNVoxel = *NZold;
  punt->xDim = *oldXps; punt->yDim = *oldYps; punt->zDim = *oldZps;
  int maxNewXVoxel,maxNewYVoxel,maxNewZVoxel;
  
  // Pulisci la matrice di destinazione
  for(int ct=0;ct< ((*pts_x) * (*pts_y) * (*pts_z)); ct++ ) returnMatrix[ct]=0;  
  ct=0;

//  for(  ct=0; ct<*numberOfPts; ct++ ) {
  for( zIndex = 0; zIndex < *pts_z; zIndex++ )   {
//    printf("\n ===> %d",zIndex);
    for( yIndex = 0; yIndex < *pts_y; yIndex++ )   {
      for( xIndex = 0; xIndex < *pts_x; xIndex++ )   {
        xPos = newXpts[xIndex];
        yPos = newYpts[yIndex];
        zPos = newZpts[zIndex];
        // acquisisci i dati relativi ai voxel della vecchia matrice i cui centroidi sono
        // ai vertici del cubo da interpolare
        _c_getCubeVertex( punt , xPos , yPos , zPos ,ct );

        // reperisci i valori di tali voxel dalla vecchia matrice
        // ed effettua l'interpolazione rispetto a tali centroidi
        valoreCalcolato = _c_TrilinearInterpolation(
          values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z1 (sample value)
          punt->xInf*(*oldXps)+(*oldXps)/2, //x0
          punt->yInf*(*oldYps)+(*oldYps)/2, //y0,
          punt->zInf*(*oldZps)+(*oldZps)/2, //z0,
          *oldXps, //dx1x0,
          *oldYps, //dy1y0,
          *oldZps, //dz1z0,
          xPos, yPos, zPos);
        
        // memorizza il risultato nella nuova struttura
        returnMatrix[ ct ] = valoreCalcolato; 
        ct++;
      }
    }
  }
  free(punt);
}
 */
/*
 * NAME: newnewtrilinearInterpolator_onGivenPoints
 * NXold,NYold,NZold : numero di voxel della matrice in ingresso
 * oldXps,oldYps,oldZps : pixelSpacing lungo le 3 direzioni della matrice in ingresso
 * newXpts,newYpts,newZpts : lista dei punti su cui andare ad interpolare
 * numberOfPts : numero di punti da interpolare
 * pts_x,pts_y,pts_z :
 */
/*
void old_newnewtrilinearInterpolator_onGivenPoints(
    int *NXold, int *NYold,int *NZold,
    double *oldXps, double *oldYps, double *oldZps,
    double *newXpts, double *newYpts, double *newZpts, int *numberOfPts,
    double *values,double *returnMatrix) {
  int ct;
  int zVoxelProgressivo,yVoxelProgressivo,xVoxelProgressivo;
  double zPos,yPos,xPos,valoreCalcolato;
  struct _c_data *punt;
  
  // Alloca un puntatore alla struttura _c_data
  punt = (struct _c_data *)calloc(1,sizeof(struct _c_data));
  if( punt == NULL ) return;
  
  // calcola il nuovo passo e copia i valori in _c_data: l'idea è quella di ridurre il passaggio parametri
  // fra funzioni ed il clone delle variabili per preservare memoria 
  //punt->newXDim = *newXps;  punt->newYDim = *newYps;  punt->newZDim = *newZps;
  punt->xNVoxel = *NXold; punt->yNVoxel = *NYold; punt->zNVoxel = *NZold;
  punt->xDim = *oldXps; punt->yDim = *oldYps; punt->zDim = *oldZps;
  int maxNewXVoxel,maxNewYVoxel,maxNewZVoxel;
  
  
  printf("\n %d,%d,%d - %lf,%lf,%lf - %lf,%lf,%lf,%d",*NXold,*NYold,*NZold,*oldXps,*oldYps,*oldZps,*newXpts,*newYpts,*newZpts,*numberOfPts);
  return;
  
  
  // Pulisci la matrice di destinazione
  for(int ct=0;ct< *numberOfPts; ct++ ) returnMatrix[ct]=0;  
  ct=0;
  
  for(  ct=0; ct<*numberOfPts; ct++ ) {
    xPos = newXpts[ct];
    yPos = newYpts[ct];
    zPos = newZpts[ct];
    // acquisisci i dati relativi ai voxel della vecchia matrice i cui centroidi sono
    // ai vertici del cubo da interpolare
    _c_getCubeVertex( punt , xPos , yPos , zPos ,ct );
    //       printf("\n (%d,%d,%d) X=%lf,%lf  Y=%lf,%lf  Z=%lf,%lf",xPos,yPos,zPos,punt->xInf,punt->xSup,punt->yInf,punt->ySup,punt->zInf,punt->zSup);
    //        return;
    printf("\n %d ",ct);
    if(ct>=10000) return;
    
    // reperisci i valori di tali voxel dalla vecchia matrice
    // ed effettua l'interpolazione rispetto a tali centroidi
    valoreCalcolato = _c_TrilinearInterpolation(
      values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z0 (sample value)
            values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z1 (sample value)
                  values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z0 (sample value)
                        values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z1 (sample value)
                              values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z0 (sample value)
                                    values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z1 (sample value)
                                          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z0 (sample value)
                                                values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z1 (sample value)
                                                      punt->xInf*(*oldXps)+(*oldXps)/2, //x0
                                                      punt->yInf*(*oldYps)+(*oldYps)/2, //y0,
                                                      punt->zInf*(*oldZps)+(*oldZps)/2, //z0,
                                                      *oldXps, //dx1x0,
                                                      *oldYps, //dy1y0,
                                                      *oldZps, //dz1z0,
                                                      xPos, yPos, zPos);
    
    // memorizza il risultato nella nuova struttura
    returnMatrix[ ct ] = valoreCalcolato;        
  }
  free(punt);
}
*/

/*
 nvert :   number of vertex
 vertx :   x coords of vertex points
 verty :   y coords of vertex points
 
 */

int isThePointInsideThePoly(int nvert, double *vertx, double *verty, double testx, double testy, int fromPosition, int toPosition)
{
  int i, j, c = 0;
  
  
  for (i = fromPosition, j = fromPosition+nvert-1; i < fromPosition+nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
         (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }  
  return c;
}
/*
 * Function for calculating the <x,y,z> coords of a voxel into the space
 * starting from rows (Ny) and columns (Nx) of the 2D Matrix of doses
 * Nx	:		columns
 * Ny 	:		rows
 * DOM  :		Dicom Orientation Matrix
 */ 
struct pointInSpace get3DPosFromNxNy(int Nx, int Ny, struct DICOM_OrientationMatrix DOM) {
  struct pointInSpace pis;
  //printf("\n a11=%lf, Nx=%d Ny=%d, a12=%lf, Sx=%lf  \n",DOM.a11,Nx, Ny,  DOM.a12, DOM.Sx);
  pis.x = DOM.a11 * Nx * DOM.XpixelSpacing + DOM.a12 * Ny * DOM.YpixelSpacing + DOM.Sx;
  pis.y = DOM.a21 * Nx * DOM.XpixelSpacing + DOM.a22 * Ny * DOM.YpixelSpacing + DOM.Sy;
  pis.z = DOM.a31 * Nx * DOM.XpixelSpacing + DOM.a32 * Ny * DOM.YpixelSpacing + DOM.Sz;
  return pis;
}
/*
 * Function to calculate the index position of an array which refers to a 3D matrix
 * x,y,z are the coors of the matrix you are interested in 
 * nx,ny,nz  are the dimensions of the matrix along the 3 axes
 *
 */
int posDecod(int x, int y, int z, int nx, int ny, int nz) {
  return( z*ny*nx+ y*nx + x   );
}
int posDecod2D(int x, int y, int nx, int ny) {
  //return(x*ny + y);
  return(nx*y + x);
}
/*
 *  Algorithm for Point In Polygon on Oblique planes
 *  calculated over a multiple series of contours each one on a different slice.
 *	Each slice is ordered according NumSlices value.
 *	Variables to be addressed to the C code are:
 *      PIPvector       vector of Points in Polygon to be arranged in an array in R
 *	totalY: 	vector of all Y coordinates of contours (one appended to the other)
 *	totalX: 	vector of all X coordinates of contours (one appended to the other)
 *      numberOfPoints: number of points declared in totalX/Y
 *      nX, nY, nZ:	x,y,z dimensions of the voxel matrix
 *      arrayAssociationROIandSlice: associates each point in totalX/Y with a slice. terminator (-100000) means the ROI is closed
 *	arrayDOM:		DICOM orientation vector, vectorized DICOM orientation matrix multiplied by pixel spacing
 * 
 *      minX,maxX,minY,maxY: min, max values of coords in order to preserve time. In future they will be calculated directly into this C code...
 */
void NewMultiPIPObl (int *PIPvector, double *totalX, double *totalY, int *numberOfPoints,
                     int *nX, int *nY, int *nZ,
                     int *arrayAssociationROIandSlice, 
                     double *arrayDOM,
                     double *minX, double *maxX, double *minY, double *maxY) {
  int x, y, z, c, fromPosition, toPosition;
  struct pointInSpace point; 
  struct DICOM_OrientationMatrix DOM;
  
  int baseCursor = 0;
  for (baseCursor = 0; baseCursor < *numberOfPoints; baseCursor++ ) {
    
    // check if a new ROI is beginning (AND you are not looking the last point!) )
    if ( totalX[ baseCursor ] == -10000 && baseCursor < (*numberOfPoints-1) ) {
      
      // get the corresponding z position
      z = arrayAssociationROIandSlice[ baseCursor + 1 ];
      
      // get the DOM
      DOM.a11= arrayDOM[ z * 9 + 0 ]; DOM.a21= arrayDOM[ z * 9 + 1 ];
      DOM.a31= arrayDOM[ z * 9 + 2 ]; DOM.a12= arrayDOM[ z * 9 + 3 ];
      DOM.a22= arrayDOM[ z * 9 + 4 ]; DOM.a32= arrayDOM[ z * 9 + 5 ];
      DOM.Sx=  arrayDOM[ z * 9 + 6 ]; DOM.Sy=  arrayDOM[ z * 9 + 7 ];
      DOM.Sz=  arrayDOM[ z * 9 + 8 ];
      DOM.XpixelSpacing = 1;  DOM.YpixelSpacing = 1;
      
      // now Run to find the beginning (clear) and the end of the 
      // points of the ROI, in the given array
      fromPosition = baseCursor;
      for ( toPosition = fromPosition+1; totalX[toPosition]!=-10000; toPosition++) ;;
      toPosition--;                
      
      // check all the voxel on such z-slice
      // loop through Y axis
      for ( y = 0; y < *nY; y++ ) {	
        // loop through X axis
        for ( x = 0; x < *nX; x++ ) {
          point = get3DPosFromNxNy(x, y, DOM);
          
          // check the box, just to avoid redundant computation
          if(point.x < *minX || point.x > *maxX || point.y < *minY || point.y > *maxY) continue;  
          
          // now check if the given point is in the poly
          c = isThePointInsideThePoly(toPosition-fromPosition-1, totalX, totalY, 
                                      point.x, point.y,  fromPosition+1,  toPosition+1);
          
          if( c == 1 )   {                                             
            PIPvector[z * (*nX) * (*nY) + x + y * (*nX)] = !PIPvector[z * (*nX) * (*nY) + x + y * (*nX)];
          }                                
        }
      }
    }
  }
} 
/*
 *  Algorithm for eroding margins in voxelcube structures
 *  cube : a pointer to the array of the cube structure
 *  nX,nY,nZ : dimensions of the cube
 *  mx,my,mz : erosion along the three axis
 *  iterator : must be set to '0', it is used to check the deep of the iteration tree.
 *  minValue : which is the minimum value that should be considered 'water' *  
 */
void erosion( double *cube, 
              int *nX, int *nY, int *nZ, 
              int *mx, int *my, int *mz, 
              int *iterator, int *minValue) {
  int x,y,z,center,ct;
  if( *iterator >= 10) return;  // just to avoid infinite loops
  if(*mx == 0 && *my ==0 && *mz ==0 ) return;
  
  
  // loop per ogni elemento del cubo
  for( z=0; z<*nZ; z++ ) {
    for( y=0; y<*nY; y++ ) {
      for( x=0; x<*nX; x++) {
        // prendi l'offset relativo al punto in esame
        center = posDecod(x,y,z,*nX,*nY,*nZ);
        // se != 0 vediamo l'intorno
        if(cube[center]>(*minValue+1)) {
          if( *mx>0 ){
            //if(x==0 || x==*nX) cube[center]=-1;
            if(x==0 || x>=(*nX-1)) cube[center]=-1;
            else if( cube[posDecod(x-1,y,z,*nX,*nY,*nZ)]<(*minValue+1) || cube[posDecod(x+1,y,z,*nX,*nY,*nZ)]<(*minValue+1)  ) cube[center]=-1;
          }
          if( *my>0 ){
            // if(y==0 || y==*nY) cube[center]=-1;
            if(y==0 || y==(*nY-1)) cube[center]=-1;
            else if( cube[posDecod(x,y-1,z,*nX,*nY,*nZ)]<(*minValue+1) || cube[posDecod(x,y+1,z,*nX,*nY,*nZ)]<(*minValue+1)  ) cube[center]=-1;
          }        
          if( *mz>0 ){
            // if(z==0 || z==*nZ) cube[center]=-1;
            if(z==0 || z==(*nZ-1)) cube[center]=-1;
            else if( cube[posDecod(x,y,z-1,*nX,*nY,*nZ)]<(*minValue+1) || cube[posDecod(x,y,z+1,*nX,*nY,*nZ)]<(*minValue+1)  )  cube[center]=-1;
          }
        }
      }
    }
  }
  // decrementa i vincoli sui margini
  if( *mx>0 ) *mx = *mx - 1;
  if( *my>0 ) *my = *my - 1;  
  if( *mz>0 ) *mz = *mz - 1;
  // trasforma tutti i '-1' in '0' per l'iterazione successiva
  //for(ct=0; ct<= posDecod((*nX-1),(*nY-1),(*nZ-1),*nX,*nY,*nZ); ct++) {
  for(ct=0; ct<= ((*nX)*(*nY)*(*nZ)-1); ct++) {
    if(cube[ct]==-1) cube[ct]=(*minValue);
  }
  // rilancia ricorsivamente
  *iterator = *iterator + 1;
  erosion( cube, nX, nY, nZ, mx, my, mz, iterator, minValue );
}

/*
 * Function to calculate the raw surface from a 3D voxel matrix
 + arr the input matrix, serialized
 * nX,nY,nZ  are the dimensions of the matrix along the 3 axes
 * pSX, pSY,pSZ are the pixelSpacing along the axec
 * surface is the output
 * NA_Val is the numerical value which represent 'NA' (the '0' in the old implementation)
 */
void rawSurface(double *arr, int *nX, int *nY, int *nZ, double *pSX, double *pSY, double *pSZ, double *surface, double *NA_Val) {
  int z,y,x;
  // reset surface value
  *surface = 0;
  *NA_Val = *NA_Val + 10;
  // loop the 3D-matrix
  for( z = 0 ; z < *nZ ; z++ ) {
    for( y = 0 ; y < *nY ; y++ ) {
      for( x = 0 ; x < *nX ; x++ ) {
        if( arr[  posDecod(x,y,z,*nX,*nY,*nZ) ] > *NA_Val ) {
          // if a non-zero voxel is on the border, surface cannot be calculated
          // for now skipped but you now... in the future....
          //        if( x==0 || x==*nX || y==0 || y==*nY || z==0 || z==*nZ) {*surface = -1; return; }
          
          // is it a border-voxel in each possible direction?
          if(z+1!=*nZ) {if( arr[  posDecod(x,y,z+1,*nX,*nY, *nZ ) ] <= *NA_Val ) *surface+= (*pSX) * (*pSY);}
          if(z-1>0) {if( arr[  posDecod(x,y,z-1,*nX,*nY, *nZ ) ] <= *NA_Val ) *surface+= (*pSX) * (*pSY);}
          if(y+1!=*nY) {if( arr[  posDecod(x,y+1,z,*nX, *nY, *nZ ) ] <= *NA_Val ) *surface+= (*pSX) * (*pSZ);}
          if(y-1>0) {if( arr[  posDecod(x,y-1,z,*nX,*nY, *nZ ) ] <= *NA_Val ) *surface+= (*pSX) * (*pSZ);}
          if(x+1!=*nX) {if( arr[  posDecod(x+1,y,z, *nX, *nY, *nZ ) ] <= *NA_Val ) *surface+= (*pSY) * (*pSZ);}
          if(x-1>0) {if( arr[  posDecod(x-1,y,z, *nX, *nY, *nZ ) ] <= *NA_Val ) *surface+= (*pSY) * (*pSZ);        }
          
          if(x-1==0) *surface+= (*pSY) * (*pSZ);
          if(x+1==*nX) *surface+= (*pSY) * (*pSZ); 
          if(y-1==0) *surface+= (*pSX) * (*pSZ);
          if(y+1==*nY) *surface+= (*pSX) * (*pSZ);          
          if(z-1==0) *surface+= (*pSX) * (*pSY);
          if(z+1==*nZ) *surface+= (*pSX) * (*pSY);         
        }
      }      
    }
  }
  return;
}
/*
 *
 */
void executeCMDLine( char **stringa, int *strLength ) {
  int i, res;
  for( i = 0; i< (*strLength) ; i++ ) {
    if( stringa[0][i] == '/' ) stringa[0][i] = '\\';
  }
  res = system( stringa[0] ) ;
}
/*
 *
 */
double somma(double *vector, int num) {
  int i=0;
  double sum;
  //Computation of total
  sum = 0;
  for (i = 0; i < num; i++) {
    sum = sum + vector[i];    
  }
  return (sum);
}
//funcCalcolaCube è la funzione che calcola il cubetto su cui andrò a calcolare la dimensione frattale
//Cioè la funzione di puntamento, che mi punta il mio cubetto nella matriciona
/*
 * le variabili in ingresso sono: 
 * arrayRM2 cioè la RM di partenza con le relative dimensioni dim1,dim2,dim3
 * CuboEntropia cioè il Cubetto che funge da kernel e spazzola RM2 di semilato ux,uy,uz  
 * MatriceEntropia cioè la matrice delle stesse dimensioni di arrayRM2 con tutti i valori di Entropia 
 * 
 */
void funcCalcolaCube(int i, int j, double *arrayRM2, double *entropyCube, int ux, int uy,int dim1, int dim2) 
{
  int s,t;
  for( t=j-uy; t<=j+uy; t++) 	{
    for( s=i-ux; s<=i+ux; s++)  {
      entropyCube[posDecod2D(s-(i-ux),t-(j-uy), (2*ux)+1, (2*uy)+1)]=arrayRM2[posDecod2D(s,t,dim1,dim2)];
    }
  }
  return;
}

void fractal2D(double *arrayRM2, int *dim1, int *dim2, double *entropyCube,int *ux, int *uy, int *N)
{
  
  //Inizializzazione delle variabili
  int i,j;
  int length; //length è la lunghezza di entropyCube
  int NumGrigi;
  int LatoCubo;
  int numeroIterazioni = 1;
  
  length=(2*(*ux)+1) * (2*(*uy)+1);
  NumGrigi=(*N);
  LatoCubo=2*(*ux);
  
  //for (i=0; i<(*dim1) * (*dim2); i=i+1) 
  //{
  // printf("\n i=%d,j=%d,arrayRM2=%lf",i,j,*(arrayRM2+i)  );
  
  //}
  
  /////for(i=*ux;i<=(*dim1 - *ux - 1);i=i+(2*(*ux)+1))
  for(j=*uy;j<=(*dim2 - *uy - 1);j=j+(2*(*uy)+1))    
  {
    /////for(j=*uy;j<=(*dim2 - *uy - 1);j=j+(2*(*uy)+1))
    for(i=*ux;i<=(*dim1 - *ux - 1);i=i+(2*(*ux)+1))
    {
      // funcCalcolaCube mi inizializza l'entropyCube ogni volta con i valori di arrayRM2 che mi servono
      //printf("\n i=%d,j=%d,arrayRM2=%lf",i,j,arrayRM2[posDecod(i,j,(*dim1),(*dim2))]);
      funcCalcolaCube( i, j, arrayRM2, entropyCube, *ux, *uy,*dim1,*dim2);
      // -im      
      //      if(numeroIterazioni==14) return;
      // -fm      
      numeroIterazioni++;
      if (somma(entropyCube,length)!=0) 
      {
        NumGrigi=NumGrigi+1;         
      }
      //printf("\n entropyCube=%lf, NumGrigi=%d", *entropyCube, NumGrigi);
    }
  }
  //printf("\n NumGrigi=%d, latocubo=%d", NumGrigi, LatoCubo);
  *N = NumGrigi;
}

/*imported from Nic */

/*
 * Function for calculating the DOUBLE of area of a facet in a triangular mesh
 */
double FacetSurface(double p1X, double p1Y, double p1Z, 
                    double p2X, double p2Y, double p2Z, double p3X, double p3Y, double p3Z) {
  double ax = p2X - p1X;
  double ay = p2Y - p1Y;
  double az = p2Z - p1Z;
  double bx = p3X - p1X;
  double by = p3Y - p1Y;
  double bz = p3Z - p1Z;
  double cx = ay*bz - az*by;
  double cy = az*bx - ax*bz;
  double cz = ax*by - ay*bx;
  //printf("\np2X*p3Y-p3X*p2Y=%lf",p2X * p3Y - p3X * p2Y );
  return sqrt(cx*cx + cy*cy + cz*cz);
}

/*
 * Function for calculating SIX TIMES the signed volume
 * of a tetrahedron in a triangular mesh
 */
double SignedVolumeOfTriangle(double p1X, double p1Y, double p1Z, 
                              double p2X, double p2Y, double p2Z, double p3X, double p3Y, double p3Z) {
  double v321 = p3X*p2Y*p1Z;
  double v231 = p2X*p3Y*p1Z;
  double v312 = p3X*p1Y*p2Z;
  double v132 = p1X*p3Y*p2Z;
  double v213 = p2X*p1Y*p3Z;
  double v123 = p1X*p2Y*p3Z;
  //printf("\ndVolume=%lf", (double)(1.0/6.0));
  return (-v321 + v231 + v312 - v132 - v213 + v123);	
}

/* 
 * Function for calculating the mesh surface
 */
void MeshSurface(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Surface) {
  int n;			// counter
  *Surface = 0;	// initial surface
  for (n=0; n<*numT; n++) {
    *Surface = *Surface + FacetSurface(X[V1[n]], Y[V1[n]], Z[V1[n]], X[V2[n]], Y[V2[n]], Z[V2[n]], X[V3[n]], Y[V3[n]], Z[V3[n]]);
  }
  *Surface=0.5 * *Surface;
}
/* 
 * Function for calulating the volume of a mesh
 */
void MeshVolume(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Volume) {
  int n;			// counter
  *Volume=0; 		// initial volume
  for (n=0; n<*numT; n++) {
    *Volume = *Volume + SignedVolumeOfTriangle(X[V1[n]], Y[V1[n]], Z[V1[n]], X[V2[n]], Y[V2[n]], Z[V2[n]], X[V3[n]], Y[V3[n]], Z[V3[n]]);
    //printf("\nX = %lf Y = %lf Z = %lf n = %d", X[V1[n]], Y[V1[n]], Z[V1[n]], n);
  }
  *Volume = fabs(*Volume * (double)(1.0/6.0));  // absolute value of a double
}

