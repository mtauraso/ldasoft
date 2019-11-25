//
//  GalacticBinaryMath.c
//
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryMath.h"


double chirpmass(double m1, double m2)
{
    return pow(m1*m2,3./5.)/pow(m1+m2,1./5.);
}


double fourier_nwip(double *a, double *b, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=0; i<n; i++)
    {
        j = i * 2;
        k = j + 1;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/Sn[i];
    }
    
    return(4.0*arg);
}

double snr(struct Source *source, struct Noise *noise)
{
    double snr2=0.0;
    switch(source->tdi->Nchannel)
    {
        case 1: //Michelson
            snr2 += fourier_nwip(source->tdi->X,source->tdi->X,noise->SnX,source->tdi->N);
            break;
        case 2: //A&E
            snr2 += fourier_nwip(source->tdi->A,source->tdi->A,noise->SnA,source->tdi->N);
            snr2 += fourier_nwip(source->tdi->E,source->tdi->E,noise->SnE,source->tdi->N);
            break;
    }
    
    return(sqrt(snr2));
}

double waveform_match(struct Source *a, struct Source *b, struct Noise *noise)
{
  int N = a->tdi->N;
  int NFFT = 2*N;
  double match=0;

  double *a_A = malloc(NFFT*sizeof(double));
  double *a_E = malloc(NFFT*sizeof(double));
  double *b_A = malloc(NFFT*sizeof(double));
  double *b_E = malloc(NFFT*sizeof(double));
  for(int i=0; i<NFFT; i++)
  {
    a_A[i] = 0.0;
    a_E[i] = 0.0;
    b_A[i] = 0.0;
    b_E[i] = 0.0;
  }

  int qmin = a->qmin - a->imin;


  //Align waveforms into arrays for summing
  for(int i=0; i<a->BW; i++)
  {
    int j = i+a->qmin-qmin;
   // printf("data->qmin=%i,i=%i, imin=%i, qmin=%i, j=%i\n",qmin,i,a->imin,a->qmin,j);
    
    if(j>-1 && j<N)
    {
      int i_re = 2*i;
      int i_im = i_re+1;
      int j_re = 2*j;
      int j_im = j_re+1;
      
      a_A[j_re] = a->tdi->A[i_re];
      a_A[j_im] = a->tdi->A[i_im];
      a_E[j_re] = a->tdi->E[i_re];
      a_E[j_im] = a->tdi->E[i_im];
    }//check that index is in range
  }//loop over waveform bins

  //Align waveforms into arrays for summing
  for(int i=0; i<b->BW; i++)
  {
    int j = i+b->qmin-qmin;
    
    if(j>-1 && j<N)
    {
      int i_re = 2*i;
      int i_im = i_re+1;
      int j_re = 2*j;
      int j_im = j_re+1;
      
      b_A[j_re] = b->tdi->A[i_re];
      b_A[j_im] = b->tdi->A[i_im];
      b_E[j_re] = b->tdi->E[i_re];
      b_E[j_im] = b->tdi->E[i_im];
    }//check that index is in range
  }//loop over waveform bins

  
  double aa = fourier_nwip(a_A,a_A,noise->SnA,N) + fourier_nwip(a_E,a_E,noise->SnE,N);
  double bb = fourier_nwip(b_A,b_A,noise->SnA,N) + fourier_nwip(b_E,b_E,noise->SnE,N);
  double ab = fourier_nwip(a_A,b_A,noise->SnA,N) + fourier_nwip(a_E,b_E,noise->SnE,N);
  
  match = ab/sqrt(aa*bb);

  free(a_A);
  free(a_E);
  free(b_A);
  free(b_E);
  
  return match;
}

// Recursive binary search function.
// Return nearest smaller neighbor of x in array[nmin,nmax] is present,
// otherwise -1
int binary_search(double *array, int nmin, int nmax, double x)
{
  int next;
  if(nmax>nmin)
  {
    int mid = nmin + (nmax - nmin) / 2;

    //find next unique element of array
    next = mid;
    while(array[mid]==array[next]) next++;
    
    // If the element is present at the middle
    // itself
    if (x > array[mid])
    {
      if(x < array[next]) return mid;
    }
    
    // the element is in the lower half
    if (array[mid] > x)
      return binary_search(array, nmin, mid, x);
    
    // the element is in upper half
    return binary_search(array, next, nmax, x);
  }
  
  // We reach here when element is not
  // present in array
  return -1;
}

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N)
{
    int i,j;
    
    // Don't let errors kill the program (yikes)
    gsl_set_error_handler_off ();
    int err=0;
    
    // Find eigenvectors and eigenvalues
    gsl_matrix *GSLfisher = gsl_matrix_alloc(N,N);
    gsl_matrix *GSLcovari = gsl_matrix_alloc(N,N);
    gsl_matrix *GSLevectr = gsl_matrix_alloc(N,N);
    gsl_vector *GSLevalue = gsl_vector_alloc(N);
    
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            if(matrix[i][j]!=matrix[i][j])fprintf(stderr,"GalacticBinaryMath.c:83: WARNING: nan matrix element, now what?\n");
            gsl_matrix_set(GSLfisher,i,j,matrix[i][j]);
        }
    }
    
    // sort and put them into evec
    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    err += gsl_eigen_symmv (GSLfisher, GSLevalue, GSLevectr, workspace);
    err += gsl_eigen_symmv_sort (GSLevalue, GSLevectr, GSL_EIGEN_SORT_ABS_ASC);
    
    // eigenvalues destroy matrix
    for(i=0; i<N; i++) for(j=0; j<N; j++) gsl_matrix_set(GSLfisher,i,j,matrix[i][j]);
    
    err += gsl_linalg_LU_decomp(GSLfisher, permutation, &i);
    err += gsl_linalg_LU_invert(GSLfisher, permutation, GSLcovari);
    
    if(err>0)
    {
        /*
         fprintf(stderr,"GalacticBinaryMath.c:98: WARNING: singluar matrix, treating matrix as diagonal\n");
         fflush(stderr);*/
        for(i=0; i<N; i++)for(j=0; j<N; j++)
        {
            evector[i][j] = 0.0;
            if(i==j)
            {
                evector[i][j]=1.0;
                evalue[i]=1./matrix[i][j];
            }
        }
        
    }
    else
    {
        
        //unpack arrays from gsl inversion
        for(i=0; i<N; i++)
        {
            evalue[i] = gsl_vector_get(GSLevalue,i);
            for(j=0; j<N; j++)
            {
                evector[i][j] = gsl_matrix_get(GSLevectr,i,j);
                if(evector[i][j] != evector[i][j]) evector[i][j] = 0.;
            }
        }
        
        //for(i=0;i<N-1;i++)for(j=i+1;j<N;j++) gsl_matrix_set(GSLcovari,j,i, gsl_matrix_get(GSLcovari,i,j) );
        
        //copy covariance matrix back into Fisher
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                matrix[i][j] = gsl_matrix_get(GSLcovari,i,j);
            }
        }
        
        //cap minimum size eigenvalues
        for(i=0; i<N; i++)
        {
            if(evalue[i] != evalue[i] || evalue[i] <= 10.0) evalue[i] = 10.;
        }
    }
    
    gsl_vector_free (GSLevalue);
    gsl_matrix_free (GSLfisher);
    gsl_matrix_free (GSLcovari);
    gsl_matrix_free (GSLevectr);
    gsl_eigen_symmv_free (workspace);
    gsl_permutation_free (permutation);
}

void invert_matrix(double **matrix, int N)
{
    int i,j;
    
    // Don't let errors kill the program (yikes)
    gsl_set_error_handler_off ();
    int err=0;
    
    // Find eigenvectors and eigenvalues
    gsl_matrix *GSLmatrix = gsl_matrix_alloc(N,N);
    gsl_matrix *GSLinvrse = gsl_matrix_alloc(N,N);
    
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            if(matrix[i][j]!=matrix[i][j])fprintf(stderr,"GalacticBinaryMath.c:172: WARNING: nan matrix element, now what?\n");
            gsl_matrix_set(GSLmatrix,i,j,matrix[i][j]);
        }
    }
    
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    err += gsl_linalg_LU_decomp(GSLmatrix, permutation, &i);
    err += gsl_linalg_LU_invert(GSLmatrix, permutation, GSLinvrse);
    
    if(err>0)
    {
        fprintf(stderr,"GalacticBinaryMath.c:184: WARNING: singluar matrix\n");
        fflush(stderr);
    }
    else
    {
        //copy covariance matrix back into Fisher
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                matrix[i][j] = gsl_matrix_get(GSLinvrse,i,j);
            }
        }
    }
    
    gsl_matrix_free (GSLmatrix);
    gsl_matrix_free (GSLinvrse);
    gsl_permutation_free (permutation);
}

void matrix_multiply(double **A, double **B, double **AB, int N)
{
  //AB = A*B
  
  int i,j,k;
  
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
    {
      AB[i][j] = 0.0;
      for(k=0; k<N; k++)
      {
        AB[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  
}


void cholesky_decomp(double **A, double **L, int N)
{
  /*
   factorize matrix A into cholesky decomposition L.
   GSL overwrites the original matrix which we want
   to preserve
   */
  int i,j;
  gsl_matrix *GSLmatrix = gsl_matrix_alloc(N,N);

  //copy covariance matrix into workspace
  for(i=0; i<N; i++) for(j=0; j<N; j++) gsl_matrix_set(GSLmatrix,i,j,A[i][j]);
  
  //make the magic happen
  gsl_linalg_cholesky_decomp1(GSLmatrix);
  
  //copy cholesky decomposition into output matrix
  for(i=0; i<N; i++) for(j=0; j<N; j++)  L[i][j] = gsl_matrix_get(GSLmatrix,i,j);
  
  //zero upper half of matrix (copy of A)
  for(i=0; i<N; i++) for(j=i+1; j<N; j++) L[i][j] = 0.0;
  
  gsl_matrix_free (GSLmatrix);
  
}


/* ********************************************************************************** */
/*                                                                                                                              */
/*                                   Fourier Tools                                    */
/*                                                                                                                              */
/* ********************************************************************************** */

void dfour1(double data[], unsigned long nn, int isign)
{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;
    
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=gsl_sf_sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=gsl_sf_sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

void drealft(double data[], unsigned long n, int isign)
{
    void dfour1(double data[], unsigned long nn, int isign);
    unsigned long i,i1,i2,i3,i4,np3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;
    
    theta=3.141592653589793/(double) (n>>1);
    if (isign == 1) {
        c2 = -0.5;
        dfour1(data,n>>1,1);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2;i<=(n>>2);i++) {
        i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    } else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        dfour1(data,n>>1,-1);
    }
}

double power_spectrum(double *data, int n)
{
    int i,j;
    double Re, Im;
    
    i = 2*n;
    j = i+1;
    Re = data[i];
    Im = data[j];
    
    return (Re*Re + Im*Im);
}

