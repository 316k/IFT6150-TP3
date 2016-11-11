/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-B.c                            */
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

#define NAME_IMG_OUT1 "photograph_original_A"
#define NAME_IMG_OUT2 "photograph_degraded_withoutNoise_A" 
#define NAME_IMG_OUT3 "photograph_restored_withoutNoise_A"  
#define NAME_IMG_OUT4 "photograph_degraded_withNoise_A" 
#define NAME_IMG_OUT5 "photograph_restored_withNoise_A"  

void AddMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,float** mat2Rin,float** mat2Iin,int lgth,int wdth);
float difference_squared(float** M1, float** M2, int length, int width);
void landweber(float**, float**, float**, float**, float**, float**, float**, float, int, int, int);
void mult_pi_fct(float** f, int length, int width);

int main(int argc,char** argv)
{
    int nb_iterations;
    int length,width;
    float var, alpha = 1;
    int size_filter; /* taille du filtre servant a ajouter du flou a l'image d'entree */
    int i, j;
    
    float** image_r; /* image d'entree */
    float** image_i; /* image d'entree */
    float** g_r;  	 /* image degradee */
    float** g_i;
    float** f; 		 /* image restoree */
    float** f_i;
    float** filter_r;
    float** filter_i;    
   
    printf("Entrez la largeur du filtre passe bas : ");
    scanf("%d",&size_filter);

    printf("\nEntrez le nombre d'itérations pour LANDWEBER: ");
    scanf("%d",&nb_iterations);
		
    /* lire l'image d'entree */	
    image_r  = LoadImagePgm(NAME_IMG_IN, &length, &width);
    image_i  = fmatrix_allocate_2d(length, width);
    
    filter_r = fmatrix_allocate_2d(length, width);
    filter_i = fmatrix_allocate_2d(length, width);

    g_r = fmatrix_allocate_2d(length, width);
    g_i = fmatrix_allocate_2d(length, width);

    f = fmatrix_allocate_2d(length, width);
    f_i = fmatrix_allocate_2d(length, width);

    
    /*Initialisation a zero de toutes les matrices */
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
          
            if((i < size_filter/2.0 || i >= length - size_filter/2.0)  &&
               (j < size_filter/2.0 || j >= width - size_filter/2.0)) {
                filter_r[i][j] = 1.0/(size_filter * size_filter);
            } else {
                filter_r[i][j] = 0.0;
            }

            filter_i[i][j] = 0.0;
            f_i[i][j] = 0.0;
        }

    /* Ajouter du flou a l'image d'entree : g = image + flou */
    FFTDD(image_r, image_i, length, width);
    FFTDD(filter_r, filter_i, length, width);
    
    MultMatrix(g_r, g_i, image_r, image_i, filter_r, filter_i, length, width);

    // Passage de g dans le domaine spatial pour la suite
    // filter doit rester dans le domaine fréquentiel
    IFFTDD(g_r, g_i, length, width);
    IFFTDD(image_r, image_i, length, width);
    
    SaveImagePgm(NAME_IMG_OUT1, g_r, length, width);
    SaveImagePgm(NAME_IMG_OUT2, g_r, length, width);
    
    /*******************************************************/
    /* restorer l'image g (elle contient du flou mais pas  */
    /* de bruit) avec LANDWEBER.                           */
	/* N'oubliez pas d'afficher le ISNR a chaque iteration */
	/*******************************************************/

    landweber(f, f_i, g_r, g_i, filter_r, filter_i,
              image_r, alpha, nb_iterations, length, width);
    
    /*Sauvegarde des images */
    
    SaveImagePgm(NAME_IMG_OUT3, f, length, width);

    printf("Entrez la variance du bruit : ");
    scanf("%f",&var);
	
    /* Ajouter du bruit a l'image floue : g = g + bruit = image + flou + bruit (add_gaussian_noise) */
    add_gaussian_noise(g_r, length, width, var);

    SaveImagePgm(NAME_IMG_OUT4, g_r, length, width);

    /*******************************************************/
	/* restorer l'image g (elle contient du flou ainsi que */
	/* du bruit) avec LANDWEBER.                           */
	/* N'oubliez pas d'afficher le ISNR a chaque iteration */
	/*******************************************************/
		
    /*Sauvegarde des images */

    landweber(f, f_i, g_r, g_i, filter_r, filter_i,
              image_r, alpha, nb_iterations, length, width);
    
    /*Sauvegarde des images */    
    SaveImagePgm(NAME_IMG_OUT5, f, length, width);

    /*Liberation memoire*/
    free_fmatrix_2d(image_r); /* image d'entree */
    free_fmatrix_2d(image_i); /* image d'entree */
    free_fmatrix_2d(g_r);  	  /* image degradee */
    free_fmatrix_2d(g_i);
    free_fmatrix_2d(f); 	  /* image restoree */
    free_fmatrix_2d(f_i);
    free_fmatrix_2d(filter_r);
    free_fmatrix_2d(filter_i);
    
    /*retour sans probleme*/ 
    printf("\n C'est fini ... \n\n\n");
    return 0; 	 
}

void AddMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
               float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
    int i,j;
 
    for(i=0;i<lgth;i++)
        for(j=0;j<wdth;j++) {
            matRout[i][j] = mat1Rin[i][j] + mat2Rin[i][j];
         
            matIout[i][j] = mat1Iin[i][j] + mat2Iin[i][j];
        }
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

void landweber(float** f, float** f_i, // Image restaurée
               float** g_r, float** g_i, // Image bruitée
               float** filter_r, float** filter_i,
               float** image_r, float alpha,
               int nb_iterations, int length, int width) {

    // Calculs intermédiaires
    float** calcul1_r = fmatrix_allocate_2d(length, width);
    float** calcul1_i = fmatrix_allocate_2d(length, width);

    float** calcul2_r = fmatrix_allocate_2d(length, width);
    float** calcul2_i = fmatrix_allocate_2d(length, width);

    float isnr;
    int i,j;
    
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            f[i][j] = g_r[i][j];
        }

    
    float num = difference_squared(image_r, g_r, length, width);
    
    for(i=0; i<nb_iterations; i++) {
        
        FFTDD(f, f_i, length, width);
        
        FFTDD(g_r, g_i, length, width);

        // Calculs
        MultMatrix(calcul1_r, calcul1_i, filter_r, filter_i, f, f_i, length, width);

        IFFTDD(calcul1_r, calcul1_i, length, width);
        IFFTDD(g_r, g_i, length, width);
        
        Mult(calcul1_r, -1, length, width);
        Mult(calcul1_i, -1, length, width);
        
        AddMatrix(calcul2_r, calcul2_i, calcul1_r, calcul1_i, g_r, g_i, length, width); 

        FFTDD(calcul2_r, calcul2_i, length, width);
        
        MultMatrix(calcul1_r, calcul1_i, calcul2_r, calcul2_i, filter_r, filter_i, length, width);

        IFFTDD(calcul1_r, calcul1_i, length, width);
        
        IFFTDD(f, f_i, length, width);

        // multiplication par alpha
        Mult(calcul1_r, alpha, length, width);
        Mult(calcul1_i, alpha, length, width);
        
        AddMatrix(f, f_i, calcul1_r, calcul1_i, f, f_i, length, width);
        
        // Calculer le ISNR

        isnr = 10 * log10(num / difference_squared(image_r, f,   length, width));

        printf("%02d - ISNR : %lf\n", i, isnr);
    }

    mult_pi_fct(f, length, width);
}

/**
 * Multiplication par une fonction porte de support [0, 255]
 */
void mult_pi_fct(float** f, int length, int width) {
    int i, j;
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            f[i][j] = fmax(0, fmin(f[i][j], 255));
        }
}
