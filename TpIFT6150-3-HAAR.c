/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  :                                            */
/* Date    :                                            */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo3.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/  
/*------------------------------------------------*/
#define NAME_IMG_IN  "photograph"

#define NAME_IMG_OUT1 "photograph_haar_C"
#define NAME_IMG_OUT2 "photograph_haar_inverse_C"

int main(int argc,char** argv)
 {
  int i,j,k,l;
  int length,width;
  int M=3; /* trois niveaux de decomposition */
  
  float** image;
  float** haar;
  float** haar_inverse;
  float** tmp1;
   
  /*Allocation memoire pour les matrices images*/
  image=LoadImagePgm(NAME_IMG_IN,&length,&width);
  haar=fmatrix_allocate_2d(length,width);
  haar_inverse=fmatrix_allocate_2d(length,width);
  tmp1=fmatrix_allocate_2d(length,width);
	
  for(i=0;i<length;i++)  
    for(j=0;j<width;j++) 
    { 
			tmp1[i][j]=0.0;
			}

	/* Transformee de Haar */
  haar2D_complete(image,haar,M,length,width);

  Recal_haar(haar,M,tmp1,length,width);
  SaveImagePgm(NAME_IMG_OUT1,tmp1,length,width);

	/* Transformee inverse de Haar */
  ihaar2D_complete(haar,haar_inverse,M,length,width);
  SaveImagePgm(NAME_IMG_OUT2,haar_inverse,length,width);

  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(image);
  free_fmatrix_2d(haar);
  free_fmatrix_2d(haar_inverse);
  free_fmatrix_2d(tmp1);
  
  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}

