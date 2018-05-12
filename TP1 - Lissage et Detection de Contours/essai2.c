#include "pgm.h"
#include <math.h>





/* av[1] contient le nom de l'image, av[2] le nom du resultat . */
int main (int ac, char **av) {
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;

  if (ac < 3) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
  /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }



  /* Calcul de la FFT de l'image tout d'abord */
  double**im3=imuchar2double(im1,nl,nc);
  oldnl=nl; oldnc=nc;
  /*  la fft demande des puissances de 2. On padde avec des 0, mais les dimensions nl et nc changent */
  im7=padimdforfft(im3,&nl,&nc);
  /*
	On peut aussi directement utiliser
	im7=padimucforfft(im1,&nl,&nc);
	sans convertir im1 en image de réels
  */

  /* Création des images pour les parties réelles et imaginaires des fft  */
  im4=alloue_image_double(nl,nc); im5=alloue_image_double(nl,nc); im6=alloue_image_double(nl,nc);

  /* mouveau code */

  /* Création de la TF du filtre gaussien G (directement) */
  double theta = 0.5;
  double pi = 3.14;
  /* Partie reelle */
  double **imGaussR = alloue_image_double(nl,nc);
  double **imGaussIm = alloue_image_double(nl,nc);
  for (int i=0; i<nl; i++) {
          for (int j=0; j<nc; j++) {
                  /* expression de la TF d'une gaussienne */
                  /* On utilise la formule discrétisée,
                  * attention aux divisions entieres
                  */
                  (*imGaussR)[i*nc + j] = exp(-2*pow(pi,2)*pow(theta,2) *
                                                    (pow((i-nc/2)/nc,2)+pow((j-nl/2)/nl,2)));
          }
  }

  /* Calcul de la fft de I (de l'image)*/
  /* fft de im7,im4 */
  fft(im7,im4,im5,im6,nl,nc);
  /* On shift l'image obtenue */
  fftshift( im5, im6, im5, im6, nl, nc);

  /* Produit de FFT(I)xFFT(G)
  * (produits de complexes)
  */
  double **imgProdR = alloue_image_double(nl,nc);
  double **imgProdIm = alloue_image_double(nl,nc);
  for (int i=0; i<nl; i++) {
          for (int j=0; j<nc; j++) {
                  (*imgProdR)[i*nc + j] = ((*imGaussR)[i*nc + j]) *
                                           ((*im3)[i*nc + j]);
          }
  }

  /*  On shift le produit */
  fftshift( imgProdR, imgProdIm, imgProdR, imgProdIm, nl, nc);


  /* On prend ensuite l'inverse avec ifft de ce produit */
  /* Creation des images pour les parties réelles et imaginaires des fft inverses */
  im9=alloue_image_double(nl,nc); im10=alloue_image_double(nl,nc);

  /* Calcul de la fft inverse de imgProdR,imgProdIm */
  ifft(imgProdR,imgProdIm,im9,im10,nl,nc);

  /* Conversion en entier8bits de la partie reelle de la fftinverse,
	  Suppresion des 0 qui ont servi a completer en utilisant la fonction crop
	 Sauvegarde au format pgm de cette image qui doit etre identique a 'linverse video
	  car on a realise la suite fftinv(fft(image))*/
  ecritureimagepgm(av[2],crop(imdouble2uchar(im9,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
  //ecritureimagepgm(av[2],imdouble2uchar(im5,nl,nc),oldnl,oldnc); marche pas
  return 0;
}
