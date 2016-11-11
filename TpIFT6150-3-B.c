/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  : Mercedes Gauthier, Nicolas Hurtubise       */
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

#define NAME_IMG_OUT1 "photograph_bruite_B"
#define NAME_IMG_OUT2 "photograph_debruite_B"

float difference_squared(float** M1, float** M2, int length, int width);
void mult_pi_fct(float** f, int length, int width);

int main(int argc,char** argv){
    int i,j,k,l;
    int length,width;
    int nbLevels;
	
	float threshold;
	float var;   
	
	printf("Entrez la variance du bruit: ");
    scanf("%f",&var);
	
	printf("Entrez le nombre de niveaux a traiter : ");
    scanf("%d",&nbLevels);
		
	printf("Entrez le seuil : ");
    scanf("%f",&threshold);

	/* ouvrir l'image d'entree */
    float** image_orig = LoadImagePgm(NAME_IMG_IN, &length, &width);
    float** image = LoadImagePgm(NAME_IMG_IN, &length, &width);
    float** haar = fmatrix_allocate_2d(length, width);
    float** tmp = fmatrix_allocate_2d(length, width);
    
	/* ajouter du bruit a l'image d'entree (add_gaussian_noise) */
	add_gaussian_noise(image, length, width, var);

    SaveImagePgm(NAME_IMG_OUT1, image, length, width);

    /* debruiter l'image en seuillant les coefficients de Haar */
    haar2D_complete(image, haar, nbLevels, length, width);

    int size_image_moyenne = length / powf(2, nbLevels);
    
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            if(i < size_image_moyenne && j < size_image_moyenne)
                continue;

            if(fabs(haar[i][j]) < threshold)
                haar[i][j] = 0;
            else
                haar[i][j] = 255 * haar[i][j]/fabs(haar[i][j]) * (fabs(haar[i][j]) - threshold) / (255 - threshold);
        }

    // Transformation inverse de Haar dans `tmp`
    ihaar2D_complete(haar,tmp,nbLevels,length,width);

  
	/* afficher l'ISBN */
    float isnr = 10 * log10(difference_squared(image_orig, image, length, width) /
                            difference_squared(image_orig, tmp, length, width));

    printf("ISNR : %lf\n", isnr);
    
 	/* sauvegarder les images */
    mult_pi_fct(tmp, length, width);
    SaveImagePgm(NAME_IMG_OUT2, tmp, length, width);

    free_fmatrix_2d(image_orig);
    free_fmatrix_2d(image);
    free_fmatrix_2d(haar);
    free_fmatrix_2d(tmp);
    
    /*retour sans probleme*/ 
    printf("\n C'est fini ... \n\n\n");
    return 0; 	 
}

float difference_squared(float** M1, float** M2, int length, int width) {
    int i,j;
    float total = 0.0;
    
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            total += SQUARE(M1[i][j] - M2[i][j]);
        }

    return total;
}

/**
 * Multiplication par une fonction porte de support [0, 255]
 */
void mult_pi_fct(float** f, int length, int width) {
    int i,j;
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            f[i][j] = fmax(0, fmin(f[i][j], 255));
        }
}
