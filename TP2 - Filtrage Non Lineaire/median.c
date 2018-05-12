#include <stdio.h>
#include <stdlib.h>
#include "pgm.h"
#include <math.h>

void libere_imd(double **im) {
        free(*im);
        free(im);
}

/*
* p : pourcentile
* t : taille bloc
*/
double estimation_bruit(double **im_input, int nlIn, int ncIn,  int t, double p)
{
        double **masque = alloue_image_double(3,3);
        masque[0][0] = 0; masque[0][1] = -1/5; masque[0][2] = 0;
        masque[1][0] = -1/5; masque[1] [1] = 1; masque[1][2] = -1/5;
        masque[2][0] = 0; masque[2][1] = -1/5; masque[2][2] = 0;
        double **im_convol = alloue_image_double(nlIn, ncIn);
        double res = 0;
        int counter = 0;
        int sum = 0;

        // convolution passe haut
        for (int x=0; x<nlIn; x++) {
                for (int y=0; y<ncIn; y++) {
                        for (int p=-1; p<=1; p++) {
                                for (int q=-1; q<=1; q++) {
                                        im_convol[x][y] +=
                                         masque[p+1][q+1]*im_input[(x+p+nlIn)%nlIn][(y+q+ncIn)%ncIn];
                                }
                        }
                }
        }
        double var_locale = 0;
        double mean = 0;
        double locale_max = 0.0;
        int *histo_varlocale = malloc(sizeof(int)*9000);
        for (int x=0; x<=nlIn; x++) {
                for (int y=0; y<=ncIn; y++) {
                        mean = 0; // local mean
                        var_locale = 0;

                        for (int i=-t; i<=t; i++) {
                                for (int j=-t; j<=t; j++) {
                                        mean += im_input[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn];
                                }
                        }
                        mean /= (2*t+1)*(2*t+1);
                        //fprintf(stderr, "mean : %f \n", mean);
                        for (int i=-t; i<=t; i++) {
                                for (int j=-t; j<=t; j++) {
                                        var_locale += pow(im_input[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn] - mean, 2);
                                }
                        }

                        var_locale /= ((double) (2*t+1)*(2*t+1));
                        var_locale = sqrt(var_locale);
                        //fprintf(stderr, "value %f\n", var_locale);
                        if (var_locale > locale_max && var_locale != locale_max)  {
                                locale_max = var_locale;
                        }
                        //fprintf(stderr, "local_max : %f, var_local : %f\n", locale_max, var_locale);
                        histo_varlocale[(int) var_locale] += 1;
                }
        }
        //fprintf(stderr, "local_max : %f", locale_max);
        sum = 0;
        counter = 0;
        //fprintf(stderr, "BEFORE : sum : %d, prod : %d", sum, ((int) (p*nlIn*ncIn)));
        while ((sum < ((int) (p*nlIn*ncIn))) && (counter<9000)) {
                sum += histo_varlocale[counter];
                //fprintf(stderr, "sum : %d\n", sum);
                counter++;
        }
        libere_imd(im_convol);
        libere_imd(masque);



        free(histo_varlocale);
        return counter*1.13;
}



/*
* ncIn, nlIn: dimensions de l'image
* n : demi taille du filtre
*/
double **filtrage_median_histo(double **im_input, int ncIn, int nlIn, int n, char *resName) {
    int k=0;
    int tmpSum = 0;
    double **im_filtre = alloue_image_double(nlIn, ncIn);
    /* pour chaque ligne de l'image  (i) */
    int *histo = malloc(256*sizeof(int));
    for (int i=0; i<nlIn; i++) {
        for (int m=0; m<256; m++) {
            histo[m] = 0;
        }
        int j0=0; // démarrage sur le premier
        for (int k=-n; k<n+1; k++) {
            for (int l=-n; l<n+1; l++) {
                /* im_input image de float mais en niveuax de gris
                * on tronque les valeurs
                */
                histo[(int) im_input[(i+k+nlIn)%nlIn][(j0+l+ncIn)%ncIn]] += 1;
            }
        }

        tmpSum =0;
        k=0;
        while (tmpSum < 2*n*(n+1)) {
            tmpSum += ((int) histo[k]);
            k++;
        }
        im_filtre[i][j0] = k;

        /* mise à jour de im_input(i,j) */
        for (int j=1; j<ncIn; j++) {
            //if (j != (n+1)) {
            /* supprime de l'histo les points de la colonne j-n-1 */
            for (int p=-n; p<n+1; p++) {
                histo[(int) im_input[(p+i+nlIn)%nlIn][(j-n-1+ncIn)%ncIn]] -= 1;
            }

            /* ajoute à l'histo les points de la colone j+n */
            for (int p=-n; p<n+1; p++) {
                histo[(int) im_input[(p+i+nlIn)%nlIn][(j+n+ncIn)%ncIn]] += 1;
            }

            tmpSum =0;
            k=0;
            while (tmpSum < 2*n*(n+1)) {
                tmpSum += ((int) histo[k]);
                k++;
            }
            im_filtre[i][j] = k;
        }
    //}
    }
    ecritureimagepgm(resName, crop(imdouble2uchar(im_filtre, nlIn,ncIn),0,0,nlIn, ncIn),nlIn,ncIn);
    free(histo);
    return im_filtre;
}

double **filtrage_bilateral(double **im_input, int ncIn, int nlIn, double sigma1, double sigma2, char *resName)
{
        double **im_res = alloue_image_double(nlIn, ncIn);
        double temp = 0.0;
        int fenMedian = 3;
        double rapport = 0.0;
        double **im_median = filtrage_median_histo(im_input, ncIn, nlIn, fenMedian, "foobar.pgm");
        for (int x=0; x<nlIn; x++) {
                for (int y=0; y<ncIn; y++) {
                        temp = 0.0;
                        rapport = 0.0;
                        for (int i=(-3)*((int) sigma1); i<=3*((int) sigma1); i++) {
                                for (int j=(-3)*((int) sigma1); j<=3*((int) sigma1); j++) {
                                        rapport = exp(-(i*i+j*j)/(2*pow(sigma1, 2))) *
                                        exp(-pow(im_median[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn]-im_median[x][y],2)/(2*pow(sigma2, 2)));
                                        temp += rapport;
                                        im_res[x][y] += rapport *
                                                        im_input[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn];
                                }
                        }
                        im_res[x][y] /= temp;
                }
        }

        libere_imd(im_median);
        ecritureimagepgm(resName, crop(imdouble2uchar(im_res, nlIn,ncIn),0,0,nlIn, ncIn),nlIn,ncIn);
        return im_res;
}

double **filtrage_adapt_recur(double **im_input, int ncIn, int nlIn, int k, char *resName)
{
    double **im_res = alloue_image_double(nlIn, ncIn);
    double **im_resNext = alloue_image_double(nlIn, ncIn);
    double **wt = alloue_image_double(nlIn, ncIn);

    /* initialisation de im_res */
    for (int x=0; x<nlIn; x++) {
        for (int y=0; y<ncIn; y++) {
            im_res[x][y] = im_input[x][y];
        }
    }
    int t = 1;
    int N = 130;
    unsigned char **imcharOrig = NULL;
    double **imOrig = NULL;
    double denom = 0;
    int nl1,nc1;
    imcharOrig  = lectureimagepgm("formes2g.pgm", &nl1, &nc1);
    imOrig = imuchar2double(imcharOrig, nl1,nc1);
    double prev_psnr = 0;
    int flag = 1;
    while (t<N) {


        /* creation de w_t */
        for (int x=0; x<nlIn; x++) {
            for (int y=0; y<ncIn; y++) {
                wt[x][y] = exp(-((pow(im_res[(x+1+nlIn)%nlIn][y]-im_res[(x-1+nlIn)%nlIn][y],2))+
                                pow(im_res[x][(y+1+ncIn)%ncIn]-im_res[x][(y-1+ncIn)%ncIn],2))
                            /(2*k*k));
            }
        }


        for (int x=0; x<nlIn; x++) {
            for (int y=0; y<ncIn; y++) {
                /* mise a jour du coeff de coord [x,y] */
                denom = 0;
                im_resNext[x][y] = 0;
                for (int i=-1; i<2; i++) {
                    for (int j=-1; j<2; j++) {
                            denom += wt[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn];
                            im_resNext[x][y] += wt[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn]*im_res[(x+i+nlIn)%nlIn][(y+j+ncIn)%ncIn];
                    }
                }
                im_resNext[x][y] /= denom;
            }
        }

        /* mise a jour de im_res */
        for (int x=0; x<nlIn; x++) {
            for (int y=0; y<ncIn; y++) {
                im_res[x][y] = im_resNext[x][y];
            }
        }
        /*if (psnr_double(im_res, imOrig, nl1, nc1)) {
                flag = 0;
        }*/
        t++;
    }
    libere_imd(wt);
    libere_imd(im_resNext);
    ecritureimagepgm(resName, crop(imdouble2uchar(im_res, nlIn,ncIn),0,0,nlIn, ncIn),nlIn,ncIn);
    //fprintf(stderr, "N = %d\n", t);
    return im_res;
}

double max(double x, double y)
{
        if (x>=y) {
                return x;
        } else {
                return y;
        }
}


/*
* t : taille des regions
* r : taille des patchs
*/
double **filtrage_NI_Mean(double **im_input, int ncIn, int nlIn, int r, int t, double sigma, double h, char *resName)
{
        double **im_res = alloue_image_double(nlIn, ncIn);
        double d_P_Q = 0;
        double w_pq = 0.0;
        double sum_w = 0.0;
        int X=0;
        int Y=0;
        for (int x=0; x<nlIn; x++) {
                for (int y=0; y<ncIn; y++) {
                        sum_w = 0;
                        /* pour tout pixel de la region centre en xy */
                        for (int i=-t; i<=t; i++) {
                                for (int j=-t; j<=t; j++) {
                                        d_P_Q = 0;
                                        X = x+i;
                                        Y = y+j;
                                        for (int k=-r; k<=r; k++) {
                                                for (int l=-r; l<=r; l++) {
                                                        d_P_Q += pow((im_input[(X+k+nlIn)%nlIn][(Y+l+ncIn)%ncIn] -
                                                                  im_input[(x+k+nlIn)%nlIn][(y+l+ncIn)%ncIn]),2);
                                                }
                                        }
                                        d_P_Q /= pow((2*r+1),2);
                                        w_pq = exp(-max(d_P_Q-2*pow(sigma,2),0)/pow(h,2));
                                        im_res[x][y] += w_pq*im_input[(X+nlIn)%nlIn][(Y+ncIn)%ncIn];
                                        sum_w += w_pq;
                                }
                        }
                        im_res[x][y] /= sum_w;
                }

        }
        ecritureimagepgm(resName, crop(imdouble2uchar(im_res, nlIn,ncIn),0,0,nlIn, ncIn),nlIn,ncIn);
        return im_res;
}


int main(int ac, char **av) {
    int nb,nl,nc, oldnl,oldnc;
    unsigned char **im2=NULL,** im1=NULL;
    double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;
    if (ac < 3) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
    unsigned char **imb, **ima; int m;
    im1=lectureimagepgm(av[1],&nl,&nc);
    oldnl = nl;
    oldnc = nc;
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
    double **im3=imuchar2double(im1,nl,nc);
    im4 = imuchar2double(im1,nl,nc);
    im5 = imuchar2double(im1,nl,nc);
    im6 = imuchar2double(im1, nl, nc);
    int nl1 = 0,nc1 = 0;
    double **imOrig = NULL;
    unsigned char **imcharOrig = NULL;
    filtrage_median_histo(im3, nc, nl, 10, av[2]);
    //filtrage_adapt_recur(im3, nc, nl, 48, av[2]);
    //filtrage_bilateral(im3, nc, nl, 5, 30, av[2]);
    //filtrage_NI_Mean(im3, nc, nl, 1, 10, 5, 0.40*5, av[2])

    /*
    double **im_filtre_median = filtrage_median_histo(im5, nc, nl, 5);
    double **im_filtre_bilateral = filtrage_bilateral(im4, nc, nl, 5, 2, 2);
    double **im_filtre_recur = filtrage_adapt_recur(im3, nc, nl, 60);
    double **im_filtre_NI_mean = filtrage_NI_Mean(im3, nc, nl, 10, 1);
    ecritureimagepgm("resFiltreRecur.pgm",imdouble2uchar(im_filtre_recur,nl,nc),nl,nc);
    ecritureimagepgm("resFiltreBilateral.pgm",imdouble2uchar(im_filtre_bilateral,nl,nc),nl,nc);
    ecritureimagepgm("resFiltreMedian.pgm",imdouble2uchar(im_filtre_median,nl,nc),nl,nc);
    ecritureimagepgm("resFiltreNiMean.pgm",imdouble2uchar(im_filtre_median,nl,nc),nl,nc);
    fprintf(stderr, "estimation bruit : %f", estimation_bruit(im3, nl, nc, 15, 0.5));
    libere_imd(im4);fprintf(stdout, "CHOIX DE SIGMA1 pour le filtre bilateral \n");
    fprintf(stdout, "for img : %s\n", av[1]);
    fprintf(stdout, "FORMAT : (sigma1, psnr) a sigma2 = 5 fixe\n");
    for (int i=1; i<=10; i+=1) {
            fprintf(stdout, "(%d, %f)\n",i, psnr_double(filtrage_bilateral(im3, nc, nl, i, 10), imOrig, oldnl, oldnc));
    }
    libere_imd(im_filtre_recur);
    libere_imd(im_filtre_NI_mean);*/

    //imcharOrig  = lectureimagepgm("formes2g.pgm", &nl1, &nc1);
    //imOrig = imuchar2double(imcharOrig, nl1,nc1);
    // Code exemple choix de la fenetre pour le median


    /*
    fprintf(stdout, "CHOIX DE TAILLE FENETRE FILTRE MEDIAN\n");
    fprintf(stdout, "for img : %s\n", av[1]);
    fprintf(stdout, "FORMAT : (taille_fen, psnr)\n");
    imcharOrig  = lectureimagepgm("formes2g.pgm", &nl1, &nc1);
    imOrig = imuchar2double(imcharOrig, nl1,nc1);
    for (int i=1; i<=20; i+=2) {
            fprintf(stdout, "(%d, %f)\n",i, psnr_double(filtrage_median_histo(im3, nc, nl, i, "resMedianBoucleFen.pgm"), imOrig, oldnl, oldnc));
    }
    fprintf(stdout, "#########\n");
    // bonne fenetre
    fprintf(stdout, "(%d, %f)\n",4, psnr_double(filtrage_median_histo(im3, nc, nl, 4, "_resMedian_fen4_2gbb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n",10, psnr_double(filtrage_median_histo(im3, nc, nl, 10, "_resMedian_fen10_2gbb90.pgm"), imOrig, oldnl, oldnc));
    // fenetre trop grande
    fprintf(stdout, "(%d, %f)\n", 30, psnr_double(filtrage_median_histo(im3, nc, nl, 30, "_resMedian_fen30_2gbb30.pgm"), imOrig, oldnl, oldnc));


    fprintf(stdout, "\n\n\n");
    fprintf(stdout, "CHOIX DE SIGMA1 pour le filtre bilateral \n");
    fprintf(stdout, "for img : %s\n", av[1]);
    fprintf(stdout, "FORMAT : (sigma1, psnr) a sigma2 = 10 fixe\n");
    for (int i=1; i<=10; i+=1) {
            fprintf(stdout, "(%d, %f)\n",i, psnr_double(filtrage_bilateral(im3, nc, nl, i, 10, "tmp.pgm"), imOrig, oldnl, oldnc));
    }

    fprintf(stdout, "CHOIX DE SIGMA2 pour le filtre bilateral \n");
    fprintf(stdout, "for img : %s\n", av[1]);
    fprintf(stdout, "FORMAT : (sigma2, psnr) a sigma1 = 5 fixe\n");
    for (int i=1; i<=150; i+=10) {
            fprintf(stdout, "(%d, %f)\n",i, psnr_double(filtrage_bilateral(im3, nc, nl, 5, i, "tmp.pgm"), imOrig, oldnl, oldnc));
    }*/
    /*
    fprintf(stdout, "(%d, %f) sigma2 a 20 fixe\n",2, psnr_double(filtrage_bilateral(im3, nc, nl, 2, 20, "influence_sig_bilat_1_20.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f) sigma2 a 20 fixe\n",5, psnr_double(filtrage_bilateral(im3, nc, nl, 5, 20, "ingluence_sig_bilat_5_20.pgm"), imOrig, oldnl, oldnc));
    */
    /*
    fprintf(stdout, "\n\n\n");
    fprintf(stdout, "CHOIX DE K pour le filtre adaptatif recursif \n");
    fprintf(stdout, "for img : %s\n", av[1]);
    fprintf(stdout, "FORMAT : (k, psnr) a sigma1 \n");
    fprintf(stdout, "estim bruit : %f", estimation_bruit(im3, nl,nc, 15, 0.5));
    for (int i=8; i<=60; i+=5) {
            fprintf(stdout, "(%d, %f)\n", i, psnr_double(filtrage_adapt_recur(im3, nc, nl, i, "resAdapt_boucle.pgm"), imOrig, oldnl, oldnc));
    }
    fprintf(stdout, "######## gen images\n");
    fprintf(stdout, "(%d, %f)\n", 19, psnr_double(filtrage_adapt_recur(im3, nc, nl, 20, "_resAdapt_k20_.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n", 19, psnr_double(filtrage_adapt_recur(im3, nc, nl, 48, "_resAdapt_k48_.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n", 19, psnr_double(filtrage_adapt_recur(im3, nc, nl, 13, "_resAdapt_k13_.pgm"), imOrig, oldnl, oldnc));


    fprintf(stdout, "\n\n\n");
    fprintf(stdout, "CHOIX DE sigma pour le filtre NI-mean avec autres parametre cf table TP \n");
    fprintf(stdout, "for img : %s\n", av[1]);
    fprintf(stdout, "FORMAT : (sigma, psnr) a sigma1 \n");
    fprintf(stdout, "estim bruit : %f", estimation_bruit(im3, nl,nc, 15, 0.5));
    */
    /*
    fprintf(stdout, "(%d, %f)\n", 5, psnr_double(filtrage_NI_Mean(im3, nc, nl, 1, 10, 3, 0.40*3, "tmp.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n", 20, psnr_double(filtrage_NI_Mean(im3, nc, nl, 2, 10, 20, 0.40*20, "tmp.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n", 35, psnr_double(filtrage_NI_Mean(im3, nc, nl, 3, 17, 35, 0.35*35, "tmp.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n", 50, psnr_double(filtrage_NI_Mean(im3, nc, nl, 4, 17, 50, 0.35*50, "tmp.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "(%d, %f)\n", 80, psnr_double(filtrage_NI_Mean(im3, nc, nl, 5, 17, 80, 0.30*80, "tmp.pgm"), imOrig, oldnl, oldnc));*/


    /*
    fprintf(stdout, "2bb90 COMPARAISON PSNR POUR LES DIFFERENTS FILTRES \n");
    fprintf(stdout, "median (%d, %f)\n",9, psnr_double(filtrage_median_histo(im3, nc, nl, 9, "=resMedian_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "bilat (%d, %f)\n",5, psnr_double(filtrage_bilateral(im3, nc, nl, 5, 30, "=resbilat_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "adapt (%d, %f)\n", 40, psnr_double(filtrage_adapt_recur(im3, nc, nl, 40, "=resAdapt_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "ni_mean (%d, %f)\n", 5, psnr_double(filtrage_NI_Mean(im3, nc, nl, 1, 10, 5, 0.40*5, "=resNImean_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "COMPARAISON TEMPS POUR LES DIFFERENTS FILTRES \n");*/

    /*
    clock_t debut,fin;
    debut=clock();
    filtrage_median_histo(im3, nc, nl, 7, "=resMedian_2gs017.pgm");
    fin = clock();
    fprintf(stdout, "filtrage median : %f\n", ((double) fin-debut)/CLOCKS_PER_SEC);

    debut=clock();
    filtrage_bilateral(im3, nc, nl, 5, 40, "=resbilat_gs017.pgm");
    fin = clock();
    fprintf(stdout, "filtrage bilat : %f\n", ((double) fin-debut)/CLOCKS_PER_SEC);


    debut=clock();
    filtrage_adapt_recur(im3, nc, nl, 60, "=resAdapt_gs017.pgm");
    fin = clock();
    fprintf(stdout, "filtrage adaptatif recursif : %f\n", ((double) fin-debut)/CLOCKS_PER_SEC);


    debut=clock();
    filtrage_NI_Mean(im3, nc, nl, 1, 10, 5, 0.40*5, "=resNImean_gs017");
    fin = clock();
    fprintf(stdout, "filtrage niMean : %f\n", ((double) fin-debut)/CLOCKS_PER_SEC);*/

    /*
    fprintf(stdout, "2bb90 COMPARAISON PSNR POUR LES DIFFERENTS FILTRES \n");
    fprintf(stdout, "median (%d, %f)\n",7, psnr_double(filtrage_median_histo(im3, nc, nl, 7, "=resMedian_2gs017.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "bilat (%d, %f)\n",5, psnr_double(filtrage_bilateral(im3, nc, nl, 5, 40, "=resbilat_gs017.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "adapt (%d, %f)\n", 60, psnr_double(filtrage_adapt_recur(im3, nc, nl, 60, "=resAdapt_gs017.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "ni_mean (%d, %f)\n", 5, psnr_double(filtrage_NI_Mean(im3, nc, nl, 1, 10, 5, 0.40*5, "=resNImean_gs017"), imOrig, oldnl, oldnc));
    fprintf(stdout, "COMPARAISON TEMPS POUR LES DIFFERENTS FILTRES \n");*/

    /*
    fprintf(stdout, "2bb90 COMPARAISON PSNR POUR LES DIFFERENTS FILTRES \n");
    fprintf(stdout, "median (%d, %f)\n",9, psnr_double(filtrage_median_histo(im3, nc, nl, 9, "GresMedian_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "bilat (%d, %f)\n",5, psnr_double(filtrage_bilateral(im3, nc, nl, 5, 30, "Gresbilat_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "adapt (%d, %f)\n", 40, psnr_double(filtrage_adapt_recur(im3, nc, nl, 40, "GresAdapt_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "ni_mean (%d, %f)\n", 5, psnr_double(filtrage_NI_Mean(im3, nc, nl, 1, 10, 5, 0.40*5, "GresNImean_2bb90.pgm"), imOrig, oldnl, oldnc));
    fprintf(stdout, "COMPARAISON TEMPS POUR LES DIFFERENTS FILTRES \n");*/
    return 0;
}
