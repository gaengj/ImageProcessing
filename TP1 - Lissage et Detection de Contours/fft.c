#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"


#define SWAP(a,b) { tempr=(a);(a)=(b);(b)=tempr; }

static int fournFFT( double *data, long nn[], int ndim, int direct) {  // Direct =1 pour la transformee, -11 poru la transformee inverse
   int idim;
   unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
   unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
   double tempi,tempr;
   double theta,wi,wpi,wpr,wr,wtemp;
   
  if (direct!=1 && direct!=-1) return -1;

   // Total number of complex values.
   for (ntot=1,idim=0;idim<ndim;idim++)
      ntot *= nn[idim];
   --data;// Due to Numerical Recipies algorithm implementation :: array[1...n];
   
   nprev=1;
   // main loop over dimensions.

   for (idim=0;idim<ndim;idim++) {
      n=nn[idim];
      nrem=ntot/(n*nprev);
      ip1=nprev << 1; // multiplied by 2.
      ip2=ip1*n;
      ip3=ip2*nrem;
      i2rev=1;
      
      for (i2=1;i2<=ip2;i2+=ip1) { // Bit reversal section of the routine.
	 if (i2 < i2rev) {
	    for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	       for (i3=i1;i3<=ip3;i3+=ip2) {
		  i3rev=i2rev+i3-i2;
		  SWAP(data[i3],data[i3rev]);
		  SWAP(data[i3+1],data[i3rev+1]);
	       }
	    }
	 }
	 ibit=ip2 >> 1; // Divided by 2.
	 while (ibit >= ip1 && i2rev > ibit) {
	    i2rev -= ibit;
	    ibit >>= 1;
	 }
	 i2rev += ibit;
      }
      ifp1=ip1;
      
      while (ifp1 < ip2) {
	 ifp2=ifp1 << 1;
	 theta=direct*(M_PI*2)/(ifp2/ip1); // isign = +1 for transform
	 wtemp=sin(0.5*theta);
	 wpr=-2.0*wtemp*wtemp;
	 wpi=sin(theta);
	 wr=1.0;
	 wi=0.0;
	 for (i3=1;i3<=ifp1;i3+=ip1) {
	    for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	       for (i2=i1;i2<=ip3;i2+=ifp2) {
		  k1=i2;
		  k2=k1+ifp1;
		  tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
		  tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
		  data[k2]=data[k1]-tempr;
		  data[k2+1]=data[k1+1]-tempi;
		  data[k1] += tempr;
		  data[k1+1] += tempi;
	       }
	    }
	    wr=(wtemp=wr)*wpr-wi*wpi+wr;
	    wi=wi*wpr+wtemp*wpi+wi;
	 }
	 ifp1=ifp2;
      }
      nprev *= n;
   }
   
   return 1; 
}

#undef SWAP

/*
 * Determines the next power of 2a from num.
 */
int nextpow2( int num )
{
    if(ispowerof2(num))
        return num;
    for(unsigned i=1; ; ++i)
        if(num >> i == 0)
            return 1 << i;
}

/**
 * Gray image.
 * Images MUST have their dimensions all powers of 2 -> reshape.
 */
static int allfft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int direct, int dimx, int dimy) {
   long nn[2];
   unsigned long n;
   double *data, *pf;
   int x,y;
   // Build a vector with couples (ims_reel[i][j],im_imag[i][j]).
   // return dimensions ->all power of 2.
   // New values are set to 0.
   nn[0] = dimx;
   nn[1] = dimy;
   n= 2;
   data=pf=(double*)calloc(2*dimx*dimy,sizeof(double));
   if (!data) { printf("Erreru allocation\n"); return -1; } 
   for (y=0; y<dimy; y++) {
      for (x=0;x<dimx;x++) {
         *(pf++) = (double)ims_reel[y][x];
	 *(pf++) = (double)ims_imag[y][x];
      }
      pf+=2*(nn[0]-x);
    }
   pf+=2*((nn[1]-y)*nn[0]);
   
   fournFFT(data,nn,n,direct);
   pf=data;
   
   // Build outputs images.
   for (y=0; y<dimy; y++) {
     for (x=0;x<dimx; x++) 
        if (direct==1) { imd_reel[y][x]= *(pf++); imd_imag[y][x]= *(pf++); }
        else { imd_reel[y][x]= *(pf++)/(nn[0]*nn[1]); imd_imag[y][x]= *(pf++)/(nn[0]*nn[1]); }
	 // Discard those points.
     pf+=2*(nn[0]-x);
   }
   pf+=2*((nn[1]-y)*nn[0]);
   free(data);
   return 1;
}

int ispowerof2(int x) { return ((x != 0) && !(x & (x - 1))); }

int fft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int dimx, int dimy) {
  if (!ispowerof2(dimx) || !ispowerof2(dimy)) { printf("Pas une puissance de 2.\n"); return -1;  }
  return allfft(ims_reel,ims_imag,imd_reel,imd_imag ,1,dimx,dimy);
}

int ifft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int dimx, int dimy) {
  if (!ispowerof2(dimx) || !ispowerof2(dimy)) { printf("Pas une puissance de 2.\n"); return -1;  }
  return allfft(ims_reel,ims_imag,imd_reel,imd_imag ,-1,dimx,dimy);
}


/* Shift les cadres de ims (reeel et imag) dans imd reeel et imaghinaire) */
void fftshift( double** imsr, double** imsi, double** imdr, double** imdi, int nl, int nc ) {
   int midx =nc >> 1;
   int midy = nl >> 1;
   int finx, finy;
   int x,y;

   // Cas impair -> +1
   finx= (nc % 2 == 1)?  midx+1 : midx;
   finy= (nl % 2 == 1)?  midy+1 : midy;

   // Cadre 1 -> Cadre 4.
   for (y=0; y < finy; y++)
      for (x=0; x < finx; x++) {
	 imdr[y][x] = imsr[y+midy][x+midx];
	 imdi[y][x] = imsi[y+midy][x+midx];
      }
   
   // Cadre 2 -> Cadre 4.
   for (y=0; y < finy; y++)
      for ( x=finx; x < nc; x++) {
	 imdr[y][x] = imsr[y+midy][x-finx];
	 imdi[y][x] = imsi[y+midy][x-finx];
      }

   
   // Cadre 3 -> cadre 2.
   for (y=finy ; y< nl; y++)
      for ( x=0; x < finx; x++) {
	 imdr[y][x] = imsr[y-finy][x+midx];
	 imdi[y][x] = imsi[y-finy][x+midx];
      }

   // Cadre 4 -> Cadre 1.
   for (y=finy; y < nl; y++)
      for ( x=finx; x < nc; x++) {
	 imdr[y][x] = imsr[y-finy][x-finx];
	 imdi[y][x] = imsi[y-finy][x-finx];
      }
 }

double** padimdforfft(double** im, int* pnl, int* pnc) {
  if (ispowerof2(*pnl) && ispowerof2(*pnc)) return im;
  else { double** res=NULL; int i,j,anl=*pnl,anc=*pnc;
    *pnl=nextpow2(*pnl); *pnc=nextpow2(*pnc); 
    if( (res=alloue_image_double(*pnl,*pnc))==NULL) return NULL;
    for (i=0; i<anl; i++) for(j=0;j<anc; j++) res[i][j]=im[i][j];
    return res;
  }
}

double** padimucforfft(unsigned char** im, int* pnl, int* pnc) {
  double** res=NULL; int i,j,anl=*pnl,anc=*pnc;
  *pnl=nextpow2(*pnl); *pnc=nextpow2(*pnc); 
  if( (res=alloue_image_double(*pnl,*pnc))==NULL) return NULL;
  for (i=0; i<anl; i++) for(j=0;j<anc; j++) res[i][j]=im[i][j];
  return res;
}
