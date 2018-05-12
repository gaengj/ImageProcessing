#include "pgm.h"

	/* 
		Exemple de code avec Entrees Sortie et transformations simples d'images
		S'utilise sous la forme  "exemple tangram.pgm res.pgm"
 	*/
main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;
  
  if (ac < 3) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
	/* Calcul de son inverse video */
  double**im3=imuchar2double(im1,nl,nc);
  oldnl=nl; oldnc=nc;
	/*  la fft demande des puissances de 2. On padde avec des 0, mais les dimensions nl et nc changent */ 
  im7=padimdforfft(im3,&nl,&nc);
  /*
	On peut aussi directement utiliser 
	im7=padimucforfft(im1,&nl,&nc);
	sans convertir im1 en image de réels
  */ 
	/* Creation des images pour les parties reelles et imagianires des fft  */
  im4=alloue_image_double(nl,nc); im5=alloue_image_double(nl,nc); im6=alloue_image_double(nl,nc);
	/* Clacul de la fft de im7,im4 */
  fft(im7,im4,im5,im6,nl,nc);
	/* Creation des images pour les parties reelles et imagianires des fft inverses */
  im9=alloue_image_double(nl,nc); im10=alloue_image_double(nl,nc); 
        /* Clacul de la fft inverse de im5,im6 */
  ifft(im5,im6,im9,im10,nl,nc);
	/* Conversion en entier8bits de la partie reelle de la fftinverse, 
	   Suppresion des 0 qui ont servi a completer en utilisant la fonction crop
	   Sauvegarde au format pgm de cette image qui doit etre identique a 'linverse video 
	   car on a realise la suite fftinv(fft(image))*/
  ecritureimagepgm(av[2],crop(imdouble2uchar(im9,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
}

