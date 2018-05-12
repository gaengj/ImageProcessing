#include "pgm.h"
        /*
                Calcul de l'inverse video d'une image. Si le parametre sortie est NULL, on cree une nouvelle image.
                On parcourt toute l'image ligen par ligne, colonne par colonne, on calcule son inverse video
                et on retourne l'image ainsi modifiee
        */
unsigned char** inverse( unsigned char** sortie,  unsigned char** entree, int nl, int nc) { int i,j;
  if (sortie==NULL) sortie=alloue_image(nl,nc);
  for(i=0; i<nl; i++)
    for(j=0; j<nc; j++)
      sortie[i][j]=255-entree[i][j];
  return sortie;
}

	/*
		Exemple de code avec Entrees Sortie et transformations simples d'images
		S'utilise sous la forme  "exemple tangram.pgm res.pgm"
 	*/
main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;

  if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
	/* Calcul de son inverse video */
  im2=inverse(NULL,im1,nl,nc);
	/* Sauvegarde dans un fichier dont le nom est passe sur la ligne de commande */
  ecritureimagepgm(av[2],im2,nl,nc);
}
