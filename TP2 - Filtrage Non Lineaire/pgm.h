#ifndef _PGM
#define _PGM 1
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#define O_BINARY 0
#define IO_LEN  (1<<30) 

#ifndef FALSE
#define FALSE (0)
#define TRUE (1)
#endif

#define MAGIC_PGM       "P5\n"
#define MAGIC_PPM       "P6\n"

unsigned char** lectureimagepgm(char* fic, int* plignes, int* pcols);
void ecritureimagepgm(char* nom, unsigned char** im, int lignes, int cols);
double** lectureimagedoubleraw(char* fic, int plignes, int pcol);
void ecritureimagedoubleraw(char*, double**,int lignes, int cols);

unsigned char ** alloue_image(int nl, int nc);
double ** alloue_image_double(int nl, int nc);
void libere_image(unsigned char** im) ;

 double** imuchar2double(unsigned char **im, int nl, int nc);
unsigned char**imdouble2uchar(double** im,int nl, int nc);

int fft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int nl, int nc);
int ifft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int nl,int nc);

void fftshift( double** imsr, double** imsi, double** imdr, double** imdi, int nl, int nc );
int nextpow2( int num );
int ispowerof2(int num);
unsigned char** crop(unsigned char **im,int oi, int oj, int fi, int fj);
double** padimdforfft(double** im, int* pnl, int* pnc); 
double** padimd(double** im, int nl, int nc, int anl, int anc); 
double** padimucforfft(unsigned char ** im, int* pnl, int* pnc); 
double eqm(unsigned char **im1, unsigned char **im2,  int nl, int nc);
double psnr(unsigned char **im1, unsigned char **im2,  int nl, int nc) ;
double psnr_double(double** r, double** i, int nl, int nc);


double** norme(double** real, double** imag, int nl, int nc);
double** phase(double** real, double** imag, int nl, int nc);
