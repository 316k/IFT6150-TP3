/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-D.c                            */
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

#define NAME_IMG_OUT1 "photograph_original_C"
#define NAME_IMG_OUT2 "photograph_degraded_C" 
#define NAME_IMG_OUT3 "photograph_restaured_C"  

void AddMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,float** mat2Rin,float** mat2Iin,int lgth,int wdth);
float difference_squared(float** M1, float** M2, int length, int width);
void landweber(float**, float**, float**, float**, float**, float**, float**, float, int, int, int);
void mult_pi_fct(float** f, int length, int width);


int main(int argc,char** argv)
{
    int nb_iterations;
    int i,j,k,l;
    int length,width;
    float var;
    int size_filter;
    int nbLevels = 3;

    float** image_r; /* image d'entree */
    float** image_i; /* image d'entree */
    float** g_r;  	 /* image degradee */
    float** g_i;
    float** f; 		 /* image restoree */
    float** f_prev;
    float** f_i;
    float** filter_r;
    float** filter_i;
    float** haar;
    float** zeros;

    
    printf("Entrez la taille du filtre passe bas : ");
    scanf("%d",&size_filter);	

    printf("Entrez la variance du bruit : ");
    scanf("%f",&var);
 
    printf("Entrez le nombre d'itérations (LANDWEBER): ");
    scanf("%d",&nb_iterations);

    /* ouvrir l'image d'entree */
    /* lire l'image d'entree */	
    image_r  = LoadImagePgm(NAME_IMG_IN, &length, &width);

    SaveImagePgm(NAME_IMG_OUT1, image_r, length, width);
    
    image_i  = fmatrix_allocate_2d(length, width);
    
    filter_r = fmatrix_allocate_2d(length, width);
    filter_i = fmatrix_allocate_2d(length, width);

    g_r = fmatrix_allocate_2d(length, width);
    g_i = fmatrix_allocate_2d(length, width);

    f = fmatrix_allocate_2d(length, width);
    f_i = fmatrix_allocate_2d(length, width);
    f_prev = fmatrix_allocate_2d(length, width);

    haar = fmatrix_allocate_2d(length, width);
    zeros = fmatrix_allocate_2d(length, width);

    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            zeros[i][j] = 0;
        }


    /* ajouter du flou et du bruit a l'image d'entree (add_gaussian_noise) */
 	    
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

    add_gaussian_noise(g_r, length, width, var);
    SaveImagePgm(NAME_IMG_OUT2, g_r, length, width);

    float isnr;
    
    k = 0;

    // Configuration initiale
    /* 1er etape : deconvolution avec 'nb_iterations' de LANDWEBER */
    landweber(f, f_i, g_r, g_i, filter_r, filter_i,
              image_r, 1, nb_iterations, length, width);
    mult_pi_fct(f, length, width);

    /* 2e etape : Filtrage dans le domaine des ondelettes */
    char continuer = 1;
    int size_image_moyenne = length / powf(2, nbLevels);

    while(continuer) {

        char filename[50];
        sprintf(filename, "iter-%02d", k);
        SaveImagePgm(filename, f, length, width);
            
        for(i=0;i<length;i++)
            for(j=0;j<width;j++) {
                f_prev[i][j] = f[i][j];
            }
        
        landweber(f, f_i, g_r, g_i, filter_r, filter_i,
                  image_r, 1, 1, length, width);
        
        haar2D_complete(f, haar, nbLevels, length, width);
        
        for(i=0;i<length;i++)
            for(j=0;j<width;j++) {
                if(i < size_image_moyenne && j < size_image_moyenne)
                    continue;

                haar[i][j] = fmax(0, SQUARE(haar[i][j]) - 3 * var) / haar[i][j];
            }

        // Transformation inverse de Haar dans `f`
        ihaar2D_complete(haar, f, nbLevels, length, width);

        mult_pi_fct(f, length, width);

        // Calculer le ISNR
        isnr = 10 * log10(difference_squared(image_r, g_r, length, width) /
                          difference_squared(image_r, f,   length, width));

        printf("%02d - ISNR : %lf\n", k, isnr);

        // Condition d'arrêt
        continuer = (difference_squared(f, f_prev, length, width) /
                     difference_squared(f_prev, zeros, length, width)) > 0.00001;
        
        printf("Condition d'arrêt : %f / %f > 10⁻⁵ = %d\n",
               difference_squared(f, f_prev, length, width),
               difference_squared(f_prev, zeros, length, width),
               continuer);
        k++;
    }

    
    /* Sauvegarde des matrices sous forme d'image pgm */
    SaveImagePgm(NAME_IMG_OUT3, f, length, width);
    
    /* Liberation memoire pour les matrices */

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

    int i,j;
    
    for(i=0;i<length;i++)
        for(j=0;j<width;j++) {
            f[i][j] = g_r[i][j];
        }

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
    }

    // mult_pi_fct(f, length, width);
}

/**
 * Multiplication par une fonction porte de support [0, 255]
 */
void mult_pi_fct(float** f, int length, int width) {
    for(int i=0;i<length;i++)
        for(int j=0;j<width;j++) {
            f[i][j] = fmax(0, fmin(f[i][j], 255));
        }
}
