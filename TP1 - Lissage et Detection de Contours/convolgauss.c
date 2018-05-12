#include "pgm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* filtrage gaussien
* resName : nom du .pmg resultat
* imDoubleInput : l'entree en tableau de double
* thetaval : valeur de sigma
* nl
* nc
*/
double **convolGauss(char * resName, double **imDoubleInput, double thetaval, int nl, int nc) {
        int oldnc = nc;
        int oldnl = nl;
        double pi = 3.14;
        double theta = thetaval;
        double** im4,** im5, ** im6, ** im7, **im8, **im9, **im10;

        /*  la fft demande des puissances de 2. On padde avec des 0, mais les dimensions nl et nc changent */
        im7=padimdforfft(imDoubleInput,&nl,&nc);


        /* Création des images pour les parties réelles et imaginaires des fft  */
        im4=alloue_image_double(nl,nc); im5=alloue_image_double(nl,nc); im6=alloue_image_double(nl,nc);

        /* Création de la TF du filtre gaussien G (directement) */
        /* Partie reelle */
        double **imGaussR = alloue_image_double(nl,nc);
        double **imGaussIm = alloue_image_double(nl,nc);
        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        /* expression de la TF d'une gaussienne */
                        /* On utilise la formule discrétisée,
                        * attention aux divisions entieres
                        */
                        (*imGaussR)[i*nl + j] = exp(-2*pow(pi,2)*pow(theta,2) *
                                                          (pow((i-nl/2.0)/((double) nl),2)+pow((j-nc/2.0)/((double) nc),2)));
                }
        }

        /* Calcul de la fft de I (de l'image)*/
        /* fft de im7,im4 */
        fft(im7,im4,im5,im6,nl,nc);


        double **im5shift = alloue_image_double(nl,nc);
        double **im6shift = alloue_image_double(nl,nc);
        /* On shift l'image obtenue */
        fftshift( im5, im6, im5shift, im6shift, nl, nc);

        /* Produit de FFT(I)xFFT(G)
        * (produits de complexes)
        */
        double **imgProdR = alloue_image_double(nl,nc);
        double **imgProdIm = alloue_image_double(nl,nc);

        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        (*imgProdR)[i*nc + j] = ((*imGaussR)[i*nc + j]) *
                                                 ((*im5shift)[i*nc + j]) -
                                                 ((*im6shift)[i*nc + j]) *
                                                 ((*imGaussIm)[i*nc + j]);
                }
        }

        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        (*imgProdIm)[i*nc + j] = ((*imGaussR)[i*nc + j]) *
                                                 ((*im6shift)[i*nc + j]) +
                                                 ((*imGaussIm)[i*nc + j]) *
                                                 ((*im5shift)[i*nc + j]);
                }
        }

        double **imgProdRshift = alloue_image_double(nl,nc);
        double **imgProdImshift = alloue_image_double(nl,nc);

        /*  On shift le produit */
        fftshift( imgProdR, imgProdIm, imgProdRshift, imgProdImshift, nl, nc);


        /* On prend ensuite l'inverse avec ifft de ce produit */
        /* Creation des images pour les parties réelles et imaginaires des fft inverses */
        im9=alloue_image_double(nl,nc); im10=alloue_image_double(nl,nc);

        /* Calcul de la fft inverse de imgProdR,imgProdIm */
        ifft(imgProdRshift,imgProdImshift,im9,im10,nl,nc);

        /* Conversion en entier8bits de la partie reelle de la fftinverse,
      	  Suppresion des 0 qui ont servi a completer en utilisant la fonction crop
      	 Sauvegarde au format pgm de cette image qui doit etre identique a 'linverse video
      	  car on a realise la suite fftinv(fft(image))*/
        ecritureimagepgm(resName,crop(imdouble2uchar(im9,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
        return im9;
}

void contour_gradient(char *resName, double **imDoubleInput, int nl, int nc) {
        double **masque_prewitt;
        double **masque_prewitt_x = alloue_image_double(3,3);
        double **masque_prewitt_y = alloue_image_double(3,3);

        /* gradient en x et y */
        double **Gx = alloue_image_double(nl,nc);
        double **Gy = alloue_image_double(nl,nc);
        double **resImg = alloue_image_double(nl,nc);
        for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                        if (j==0) {
                                masque_prewitt_x[i][j] = -1;
                        } else if (j == 1) {
                                masque_prewitt_x[i][j] = 0.0;
                        } else if (j == 2) {
                                masque_prewitt_x[i][j] = 1;
                        }
                }
        }

        for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                        if (i==0) {
                                masque_prewitt_y[i][j] = -1;
                        } else if (i == 1) {
                                masque_prewitt_y[i][j] = 0.0;
                        } else if (i == 2) {
                                masque_prewitt_y[i][j] = 1;
                        }
                }
        }

        /*  M1*I  = Gx  ;  M2*I = Gy */

        /* calcul de Gx */
        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        for (int k=-1; k<=1; k++) {
                                for (int l=-1; l<=1; l++) {
                                        Gx[i][j] += masque_prewitt_x[k+1][l+1]*imDoubleInput[(i+k+nl)%nl][(j+l+nc)%nc];
                                }
                        }
                }
        }

        /* calcul de Gy */
        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        for (int k=-1; k<=1; k++) {
                                for (int l=-1; l<=1; l++) {
                                        Gy[i][j] += masque_prewitt_y[k+1][l+1]*imDoubleInput[(i+k+nl)%nl][(j+l+nc)%nc];
                                }
                        }
                }
        }

        /* Par seuil, methode la plus simple */
        double seuil_gradient = 100;
        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        if (sqrt(pow(Gx[i][j],2)+pow(Gy[i][j],2)) > seuil_gradient) {
                                resImg[i][j] = imDoubleInput[i][j];
                        }
                }
        }

        ecritureimagepgm(resName,crop(imdouble2uchar(resImg,nl,nc),0,0,nl,nc),nl,nc);
}

void contour_laplacien(char *resName, double **imDoubleInput, double theta, int nl, int nc) {
        double **fftLoGIm = alloue_image_double(nl,nc);
        double **fftLoGR = alloue_image_double(nl,nc);
        double pi = 3.14;
        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        fftLoGR[i][j] = -4*pi*(pow((i-nl/2.0)/nl,2)+pow((j-nc/2.0)/nc,2))*
                                    exp(-2*pow(pi,2)*pow(theta,2)*(pow((i-nl/2.0)/nl,2)+pow((j-nc/2.0)/nc,2)));
                }
        }

        double **imfftR = alloue_image_double(nl,nc);
        double **imfftIm = alloue_image_double(nl,nc);

        double **InputIm = alloue_image_double(nl,nc);

        /* Calcul de la fft de I (de l'image)*/
        fft(imDoubleInput,InputIm,imfftR,imfftIm,nl,nc);


        double **imfftRshift = alloue_image_double(nl,nc);
        double **imfftImshift = alloue_image_double(nl,nc);


        /* On shift l'image obtenue */
        fftshift( imfftR, imfftIm, imfftRshift, imfftImshift, nl, nc);

        /* Produit des ffts */
        double **imgProdR = alloue_image_double(nl,nc);
        double **imgProdIm = alloue_image_double(nl,nc);

        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        (*imgProdR)[i*nl + j] = ((*fftLoGR)[i*nl + j]) *
                                                 ((*imfftRshift)[i*nl + j]) -
                                                 ((*imfftImshift)[i*nl + j]) *
                                                 ((*fftLoGIm)[i*nl + j]);
                }
        }

        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        (*imgProdIm)[i*nc + j] = ((*fftLoGR)[i*nl + j]) *
                                                 ((*imfftImshift)[i*nl + j]) +
                                                 ((*fftLoGIm)[i*nl + j]) *
                                                 ((*imfftRshift)[i*nl + j]);
                }
        }

        double **imgProdRshift = alloue_image_double(nl,nc);
        double **imgProdImshift = alloue_image_double(nl,nc);

        /*  On shift le produit */
        fftshift( imgProdR, imgProdIm, imgProdRshift, imgProdImshift, nl, nc);

        /* Creation des images pour les parties réelles et imaginaires des fft inverses */
        double **laplacienR=alloue_image_double(nl,nc);
        double **laplacienIm=alloue_image_double(nl,nc);

        /* Calcul de la fft inverse de imgProdR,imgProdIm */
        ifft(imgProdRshift,imgProdImshift,laplacienR,laplacienIm,nl,nc);

        double **res = alloue_image_double(nl,nc);

        /* Maintenant on ne garde que certains points : la ou le laplacien s'annule */
        for (int i=0; i<nl; i++) {
                for (int j=0; j<nc; j++) {
                        if (((j+1) < nc) && ((i+1) < nl)) {
                                if (!((laplacienR[i][j]*laplacienR[i+1][j] > 0) &&
                                 ((laplacienR[i][j]*laplacienR[i][j+1])>0))) {
                                        res[i][j] = imDoubleInput[i][j];
                                 }
                        }
                }
        }

        ecritureimagepgm(resName,crop(imdouble2uchar(res,nl,nc),0,0,nl,nc),nl,nc);

}

/* filtrage spatial
* resName : nom du .pmg resultat
* imDoubleInput : l'entree en tableau de double
* thetaval : valeur de sigma
* nl
* nc
* W_sigma : valuer de W(sigma) choisie
* ecrit l'image ds resName
* renvoit le tableau de double resultat
*/
double **convolSpatiale(char *resName, double **imDoubleInput, double theta,int nl, int nc, int W_sigma) {
        double **im1, **im2, **im3, **im4;
        double **imR_reel = alloue_image_double(nl,nc);
        double **imR_img = alloue_image_double(nl,nc);
        double **coeff_convol_R = alloue_image_double(nl,nc);
        double **coeff_convol_Rtmp = alloue_image_double(nl,nc);
        double **coeff_convol_img = alloue_image_double(nl,nc);
        int W_theta = W_sigma;//3*((int) ceil(theta))+1; /* taille du filtre */

        /* Convolution spatiale */
        /* On precalcul les coefficients */
        double coeff_separable[nl];
        double tmpSum = 0.0;

        /* Le filtre gaussien est linéairement séparable */
        for (int i=1; i<=W_theta; i++) {
                tmpSum += exp(-pow(i,2)/((double) 2*pow(theta, 2)));
        }

        tmpSum *= 2;
        tmpSum += exp(-pow(0,2)/((double) 2*pow(theta, 2)));

        for(int i=0; i<=W_theta; i++) {
                /* 2 tmpSun car la somme va de -W a +W */
                coeff_separable[i] = exp(-pow(i,2)/((double) 2*pow(theta,2)))/tmpSum;
        }

        /* Calcul de la convolution spatiale proprement dite */
        for (int x=0; x<nl; x++) {
                for (int y=0; y<nc; y++) {
                        // convolution selon x
                        //coeff_convol_R[x][y] = 0.0;
                        for (int i=(-W_theta); i<=W_theta; i++) {
                                //printf("coeff : %f ", coeff_separable[(int) abs(i)])
                                coeff_convol_Rtmp[x][y] += coeff_separable[(int) abs(i)] *
                                        imDoubleInput[(x+i+nl)%nl][(y+nc)%nc];
                        }
                }
        }

        for (int x=0; x<nl; x++) {
                for (int y=0; y<nc; y++) {
                        // convolution selon y
                        for (int j=(-W_theta); j<=W_theta; j++) {
                                coeff_convol_R[x][y] += coeff_separable[(int) abs(j)] *
                                        coeff_convol_Rtmp[(x+nl)%nl][(y+j+nc)%nc];
                        }
                }
        }
        /*
        for (int x=0; x<nl; x++) {
                for (int y=0; y<nc; y++) {
                        coeff_convol_R[x][y] /= 1/(2*3.14*pow(theta,2));
                }
        }*/
        unsigned char **imsource = NULL;
        double **imsourcedouble = NULL;
        int mnl, mnc;
        imsource = lectureimagepgm("formes2g.pgm", &mnl, &mnc);
        imsourcedouble = imuchar2double(imsource,mnl,mnc);

        ecritureimagepgm(resName,imdouble2uchar(coeff_convol_R,nl,nc),nl,nc);
        return coeff_convol_R;
}


/* av[1] contient le nom de l'image, av[2] le nom du resultat . */
int main (int ac, char **av) {
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;

  double theta = 3;
  double pi = 3.14;



  if (ac < 3) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
  // Lecture d'une image pgm dont le nom est passe sur la ligne de commande
  printf("reading image %s\n", av[1]);
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
  double **im3=imuchar2double(im1,nl,nc);

  // exemple d'utilisation des fonctions avec des arguments ligne de commande
  // nous avons enlevé du main les boucles générants les psnr
  // et tout les tests sur des images avec différents sigma pour des raisons
  // de lisibilité
  convolGauss(av[2], im3, theta, nl, nc);
  //int w_sigma = 3*ceil(theta)+1;
  //convolSpatiale(av[2],im3,theta,nl,nc,w_sigma);

  //contour_gradient(av[2], im3, nl, nc);
  //contour_laplacien(av[2], im3, 4, nl, nc);

  // on utilise le psnr comme suit :


  return 0;
}
