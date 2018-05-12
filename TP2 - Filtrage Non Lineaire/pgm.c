#include "pgm.h"

static int match_key(int fd, char *key) { char buf[80];
  read(fd, buf, strlen(key));
  if( strncmp(buf, key, strlen(key)) != 0 ) return FALSE;
  else return TRUE;
}

static void skip_comment(int fd, char code, char *c) {
  while( *c == code ) {
      while( (read(fd, c, 1) == 1 ) && (*c != '\n') ) ;
      read(fd, c, 1);
    }
}

static void read_header_line(int fd, char *buf) { int i;
  i = 1;
  while( (read(fd, &buf[i], 1) == 1 ) && (buf[i] != '\n') && (buf[i] != '\r') && (i<79) ) i++;
  buf[i] = 0;
}

static int get_pgm_header(int fd, char *magic, int *width, int *height) { char buf[80]; 
  if( !match_key(fd, magic) ) return FALSE;
  read(fd, buf, 1);
  skip_comment(fd, '#', buf);
  read_header_line(fd, buf);
  sscanf(buf, "%d %d", width, height);
  read(fd, buf, 1);
  skip_comment(fd, '#', buf);
  read_header_line(fd, buf);
  return TRUE;
}

static int open_read(char *filename) { int fd; 
  if( (fd = open(filename, O_BINARY|O_RDONLY)) < 0 )
    fprintf(stderr, "can't reset file `%s': %s\n", filename, strerror(errno));
  return fd;
}

static int open_read_pixmap(char *filename, char *magic, int *width, int *height) { int fd; 
  if( (fd = open_read(filename)) < 0) return fd;
  if( !get_pgm_header(fd, magic, width, height) ) {
      fprintf(stderr, "can't read header of %s\n", filename);
      return -1;
    }
  return fd;
}

static unsigned char *alloc_pixmap(long size) { unsigned char *data; 
  if( (data = (unsigned char *)malloc(size)) == NULL ) {
    fprintf(stderr, "malloc error\n");
    return NULL;
    }
  return data;
}

static void load_data(int fd, unsigned char *data, long size) { char *buffer;
  int count;

  buffer = (char *)data;
  while( size > 0 ) {
    count = IO_LEN;
    if( count > size ) count = size;
    read(fd, buffer, count);
    buffer += count;
    size -= count;
  }
}

static unsigned char *load_pixmap(char *filename, int *width, int *height) { int fd;
  long size;
  unsigned char *data;

  if( (fd = open_read_pixmap(filename, MAGIC_PGM, width, height)) < 0) return NULL;
  size = (long)*width * *height;
  data = alloc_pixmap(size);
  if( data != NULL ) load_data(fd, data, size);
  close(fd);
  return data;
}

static void put_header_line(int fd, char *buf) {
  write(fd, buf, strlen(buf));
}

static void put_header_info(int fd, char *mark, char *filename) { }

static void put_pgm_header(int fd, char *magic, int width, int height, char *filename) { char buf[80]; 
  put_header_line(fd, magic);
  put_header_info(fd, "# ", filename);
  sprintf(buf, "%d %d\n255\n", width, height);
  put_header_line(fd, buf);
}

static int open_write(char *filename) { int fd; 
  if( (fd = open(filename, O_TRUNC|O_CREAT|O_BINARY|O_RDWR,0640)) < 0 )
/*
  if( (fd = open(filename, O_BINARY|O_CREAT|O_TRUNC|O_RDWR, S_IREAD|S_IWRITE)) < 0 )
*/
    fprintf(stderr, "can't Rewrite file `%s': %s\n", filename, strerror(errno));
  return fd;
}

static void store_data(int fd, unsigned char *data, long size) { char *buffer;
  int count;
  
  buffer = (char *)data;
  while( size > 0 ) {
    count = IO_LEN;
    if( count > size ) count = size;
    write(fd, buffer, count);
    buffer += count;
    size -= count;
    }
}

static void store_pixmap(char *filename, unsigned char *data, int width, int height) { int fd; 
  if( (fd = open_write(filename)) < 0 ) return;
  put_pgm_header(fd, MAGIC_PGM, width, height, filename);
  store_data(fd, data, (long)width*height);
  close(fd);
}
        /*
                Lecture dnas le fichier "filename" de format pgm d'une image en niveau de gris  de rows lignes et cols colonnes
      		Le nombre de ligne et de colonnes est lu dans le fichier "filename"
		 L'image est cree dynamiquement, les pixels remplis puis l'image est retournee
                Retourne NULL en cas d'echec
        */
unsigned char ** lectureimagepgm(char *name, int* rows, int* cols) { unsigned char ** im;
  unsigned char * p;
  int i;
  p = load_pixmap(name,cols,rows);

  if (p==NULL) return NULL;
  if ( (im=calloc(*rows,sizeof(*im))) == NULL) return NULL;
  for (*im=p,i=1; i< *rows; i++) im[i] = im[i-1]+ *cols;
  return im; 
}
        /*
                Ecriture dnas le fichier "filename" de format pgm d'une image en niveau de gris  de rows lignes et cols colonnes
        */
void ecritureimagepgm(char *name, unsigned char **im, int rows, int cols) {
  store_pixmap(name,*im,cols,rows);
}

        /*
                Ecriture dnas le fichier "filename"  d'une image de reels double precision
                Le format du fichier est "RAW", ie les donnees sont ecrites en binaire ligne par ligne et colonnes par colonnes.
        */
void ecritureimagedoubleraw(char *filename, double **im, int rows, int cols) { int fd;
  if( (fd = open_write(filename)) < 0 ) return;
  write(fd,*im,rows*cols*sizeof(**im));
  close(fd);
}

	/*
		Lecture et creation d'uen image de reels double precision contenu dns le fichier "filename"
		Le format du fichier est "RAW", ie les donnees sont ecrites en binaire ligne par ligne et colonnes par colonnes.
		Le nombre de ligne et de colonnes doivent etre connu
		L'image est cree dynamiquement, les pixels remplis
		et l'image est retournee
		Retourne NULL en cas d'echec
	*/
double** lectureimagedoubleraw( char *filename, int rows, int cols) { int fd;
  double** im=NULL;
  if( (fd = open_read(filename)) < 0 ) return NULL;
  if ( (im=alloue_image_double(rows,cols))==NULL) return NULL;
  read(fd,*im,rows*cols*sizeof(**im));
  close(fd);
}

        /*  
                Creation d'une image de d'entiers 8bits non signes de nl lignes et nl colonnes
                Rq : l'allocation des donnees est faite en une seule allocation
                        on peut donc utiliser im[i][j] pour acceder a l'element d'indice i,j 
                        ou bien *((*im)+i*nc+j) ou (*im)[i*nc+j] comme un tableau 1D
        */

unsigned char ** alloue_image(int nl, int nc) { int i;
  unsigned char** p;
  if ( (p=(unsigned char**)calloc(nl,sizeof(*p)))==NULL) return NULL;
  if ( (*p=(unsigned char*)calloc(nl*nc,sizeof(**p)))==NULL) return NULL;
  for (i=1; i<nl; i++) p[i]=p[i-1]+nc;
  return p;
}

	/*
		Creation d'une image de reelle double precision de nl lignes et nl colonnes
		Rq : l'allocation des donnees est faite en une seule allocation
			on peut donc utiliser im[i][j] pour acceder a l'element d'indice i,j 
			ou bien *((*im)+i*nc+j) ou (*im)[i*nc+j] comme un tableau 1D
	*/
double ** alloue_image_double(int nl, int nc) { int i;
  double** p;
  if ( (p=(double **)calloc(nl,sizeof(*p)))==NULL) return NULL;
  if ( (*p=(double*)calloc(nl*nc,sizeof(**p)))==NULL) return NULL;
  for (i=1; i<nl; i++) p[i]=p[i-1]+nc;
  return p;
}

	/* 	
		Libere la memoire associe a l'image im
	*/
void libere_image(unsigned char** im) {
  free(*im); free(im);
}

	/*
		Conversion d'une image de 8bits non signes en une image de double
		Utilise en particulier pour fft
	*/
double** imuchar2double(unsigned char **im, int nl, int nc) { int i;
  double** res;
  if ( (res=alloue_image_double(nl,nc))==NULL) return NULL;
  for (i=0; i<nl*nc; i++) (*res)[i]=(double) (*im)[i];
  return res;
}

	/*
		Conversion d'une image de reels double en une image 8bits non signes
	*/
unsigned char**imdouble2uchar(double** im,int nl, int nc) { int i;
   unsigned char** res;
   if ( (res=alloue_image(nl,nc))==NULL) return NULL;
   for (i=0; i<nl*nc; i++) (*res)[i]=(unsigned char) ((*im)[i]+0.5);
   return res;
}

        /*  
                Conversion d'une image de reels double en une image 8bits signes
        */
char**imdouble2char(double** im,int nl, int nc) { int i;
   char** res;
   if ( (res=alloue_image(nl,nc))==NULL) return NULL;
   for (i=0; i<nl*nc; i++) (*res)[i]=(char) ((*im)[i]+0.5);
   return res;
}

	/* 
		Copie dasn une nouvelle image de la partie de l'image comprise entre les indices (oi,oj) et (fi,fj) de l'image im 
        */
unsigned char** crop(unsigned char **im,int oi, int oj, int fi, int fj) { int i,j,nl,nc;
  unsigned char ** res;
  nl=fi-oi; nc=fj-oj;
  if ( (res=alloue_image(nl,nc))==NULL) return NULL;
  for(i=0; i<nl; i++)
     for(j=0; j<nc; j++)
       res[i][j]=im[oi+i][oj+j];
  return res;
}

