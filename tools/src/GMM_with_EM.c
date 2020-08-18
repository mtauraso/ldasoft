/*
 *  Author: Tyson B. Littenberg (MSFC-ST12)
 *  Created: 07.27.2020
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 *
 *  Library for Gaussian Mixture Model using the Expectation-Maximization Algorithm
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "GMM_with_EM.h"

static void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

void printUsage(const char *program)
{
    fprintf( stdout,"\n");
    fprintf( stdout, "Gaussian Mixture Model (GMM):\n");
    fprintf( stdout, "  Expectation-Maximization algorithm to fit GMM to \n");
    fprintf( stdout, "  input list of data samples e.g. from MCMC\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "Usage: %s required [optional] \n", program );
    fprintf( stdout, "\n");
    fprintf( stdout, "  Required:\n");
    fprintf( stdout, "     -f, --file=FILE        filename for input chain\n" );
    fprintf( stdout, "     -n, --nparams=INT      number of model parameters\n" );
    fprintf( stdout, "                            extra columns are ignored\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "  Optional:\n");
    fprintf( stdout, "    [-h, --help]            print this message and exit\n" );
    fprintf( stdout, "    [-l, --log=INT]         use log(params) in specified column number [0,m-1]\n" );
    fprintf( stdout, "                            multiple arguments can be given.\n");
    fprintf( stdout, "    [-m, --modes=INT]       number of modes in fit (2)\n");
    fprintf( stdout, "    [-s, --seed=LONG]       RNG seed\n");
    fprintf( stdout, "    [-t, --thin=INT]        downsample rate for chains (1)\n");
    fprintf( stdout,"\n");
}

void alloc_MVG(struct MVG *mode, size_t N)
{
    mode->size = N;
    mode->mu = gsl_vector_alloc(mode->size);
    mode->C = gsl_matrix_calloc(mode->size,mode->size);
    mode->L = gsl_matrix_calloc(mode->size,mode->size);
    mode->Cinv = gsl_matrix_calloc(mode->size,mode->size);
    mode->evectors = gsl_matrix_calloc(mode->size,mode->size);
    mode->evalues = gsl_vector_calloc(mode->size);
}

void free_MVG(struct MVG *mode)
{
    gsl_vector_free(mode->mu);
    gsl_matrix_free(mode->C);
    gsl_matrix_free(mode->L);
    gsl_matrix_free(mode->Cinv);
    gsl_matrix_free(mode->evectors);
    gsl_vector_free(mode->evalues);
    free(mode);
}

void write_MVG(struct MVG *mode, FILE *fptr)
{
    //pack detC,p,and Neff into a vector to make life easier
    gsl_vector *temp = gsl_vector_alloc(3);
    gsl_vector_set(temp,0,mode->detC);
    gsl_vector_set(temp,1,mode->p);
    gsl_vector_set(temp,2,mode->Neff);

    //write 'em!
    gsl_vector_fwrite(fptr,mode->mu);
    gsl_matrix_fwrite(fptr,mode->C);
    gsl_matrix_fwrite(fptr,mode->L);
    gsl_matrix_fwrite(fptr,mode->Cinv);
    gsl_matrix_fwrite(fptr,mode->evectors);
    gsl_vector_fwrite(fptr,mode->evalues);
    gsl_vector_fwrite(fptr,temp);

    gsl_vector_free(temp);
}

void read_MVG(struct MVG *mode, FILE *fptr)
{
    //vector for holding packed detC,p,and Neff
    gsl_vector *temp = gsl_vector_alloc(3);

    //read 'em!
    gsl_vector_fread(fptr,mode->mu);
    gsl_matrix_fread(fptr,mode->C);
    gsl_matrix_fread(fptr,mode->L);
    gsl_matrix_fread(fptr,mode->Cinv);
    gsl_matrix_fread(fptr,mode->evectors);
    gsl_vector_fread(fptr,mode->evalues);
    gsl_vector_fread(fptr,temp);

    //unpack 'em!
    mode->detC = gsl_vector_get(temp,0);
    mode->p = gsl_vector_get(temp,1);
    mode->Neff = gsl_vector_get(temp,2);

    gsl_vector_free(temp);
}

double multivariate_gaussian(gsl_vector *x, struct MVG *mvg)
{
    
    size_t N = x->size;
    gsl_vector *dx = gsl_vector_alloc(N);
    gsl_vector *mu = mvg->mu;
    gsl_matrix *Cinv = mvg->Cinv;
    double detC = mvg->detC;
    
    // x-mu
    for(size_t n=0; n<N; n++)
    {
        double xi = gsl_vector_get(x,n);
        double mi = gsl_vector_get(mu,n);
        double d  = xi-mi;
        
        gsl_vector_set(dx,n,d);
    }
    
    
    /*
     (x-mu)^T C^-1 (x-mu)
     */
    double CdotX;
    double chi2 = 0.0;
    for(size_t m=0; m<N; m++)
    {
        CdotX = 0.0;
        for(size_t n=0; n<N; n++)
        {
            CdotX += gsl_matrix_get(Cinv,m,n)*gsl_vector_get(dx,n);
        }
        chi2 += CdotX*gsl_vector_get(dx,m);
    }
    gsl_vector_free(dx);
    
    double p = exp(-0.5*chi2)/sqrt(pow(2.0*M_PI,N)*detC);
    
    return p > PMIN ? p : PMIN;
}

void invert_gsl_matrix(gsl_matrix *A, gsl_matrix *Ainv, gsl_matrix *L, double *detA, double *R)
{
    gsl_set_error_handler_off();
    
    //error catchers
    int err = 0;
    
    //get size of matrix (assumed to be NxN)
    size_t N = A->size1;
    
    //some workspace
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    //copy A into Ainv because LU decomposition destroys the matrix
    gsl_matrix_memcpy(Ainv,A);
    
    //cholesky decomposition
    int i;
    err += gsl_linalg_cholesky_decomp(Ainv);
    
    //get condition number
    gsl_vector *work = gsl_vector_alloc(3*N);
    err += gsl_linalg_cholesky_rcond(Ainv, R, work);

    //inverse of A
    err += gsl_linalg_cholesky_invert(Ainv);
    
    //get deteriminant, need LU decomposition
    gsl_matrix_memcpy(L,A);
    gsl_linalg_LU_decomp(L,permutation,&i);
    *detA = gsl_linalg_LU_det(L,i);
    
    //zero upper half of matrix (copy of A)
    int j;
    for(i=0; i<N; i++) for(j=i+1; j<N; j++) gsl_matrix_set(L,i,j,0.0);
    
    //clean up
    gsl_vector_free (work);
    gsl_permutation_free (permutation);
}

void decompose_matrix(gsl_matrix *A, gsl_matrix *evec, gsl_vector *eval)
{
    //error catchers
    int err = 0;
    
    //get size of matrix (assumed to be NxN)
    size_t N = A->size1;
    
    //get deteriminant, need LU decomposition
    gsl_matrix *Atemp = gsl_matrix_calloc(N,N);
    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    //copy A into Atemp because eigen_symmv destroys the matrix
    gsl_matrix_memcpy(Atemp,A);
    
    //the reason we're all here...
    err += gsl_eigen_symmv (Atemp, eval, evec, workspace);
    
    gsl_matrix_free (Atemp);
    gsl_eigen_symmv_free (workspace);
    gsl_permutation_free (permutation);
    
}

double log_likelihood(struct MVG **modes, struct Sample **samples, int NMCMC, int NMODE)
{
    
    double logL = 0.0;
    for(size_t i=0; i<NMCMC; i++)
    {
        double P = 0.0;
        for(size_t k=0; k<NMODE; k++)
        {
            P += modes[k]->p*gsl_vector_get(samples[i]->p,k);
        }
        if(P==0) exit(1);
        logL += log(P);
    }
    return logL;
}

void print_1D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMCMC, char root[], size_t ix)
{
    char filename[128];
    sprintf(filename,"%s_%i.dat",root,(int)ix);
    FILE *fptr = fopen(filename,"w");
    
    size_t NMODE = samples[0]->p->size;
    
    double *xvec = malloc(NMCMC*sizeof(double));
    double xmin,xmax;
    double x0 =  1e60;
    double xf = -1e60;
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t i=0; i<NMCMC; i++) xvec[i] = gsl_vector_get(samples[i]->x,ix);
        gsl_stats_minmax(&xmin,&xmax,xvec,1,NMCMC);
        if(xmin<x0) x0 = xmin;
        if(xmax>xf) xf = xmax;
    }
    
    double p;
    double x;
    double dx = (xf-x0)/100.;
    for(int n=0; n<=100; n++)
    {
        p = 0.0;
        x = xmin + (double)n*dx;
        for(size_t k=0; k<NMODE; k++)
        {
            double mean = gsl_vector_get(modes[k]->mu,ix);
            double var  = gsl_matrix_get(modes[k]->C,ix,ix);
            p += modes[k]->p*exp( -0.5*(x-mean)*(x-mean)/var )/sqrt(2*M_PI*var);
        }
        fprintf(fptr,"%.16g %.16g\n",x,p);
    }
    
    fclose(fptr);
    
    
    free(xvec);
}

void print_2D_contours(struct MVG **modes, size_t NMODE, char root[], size_t x1, size_t x2)
{
    char filename[128];
    FILE *fptr = NULL;
    
    struct MVG **submodes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t k=0; k<NMODE; k++)
    {
        submodes[k] = malloc(sizeof(struct MVG));
        alloc_MVG(submodes[k], 2);
    }
    
    //Pick parameters
    size_t X[2] = {x1,x2};
    
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t n=0; n<2; n++)
        {
            gsl_vector_set(submodes[k]->mu,n,gsl_vector_get(modes[k]->mu,X[n]));
            for(size_t m=0; m<2; m++) gsl_matrix_set(submodes[k]->C,m,n,gsl_matrix_get(modes[k]->C,X[m],X[n]));
        }
        
        decompose_matrix(submodes[k]->C, submodes[k]->evectors, submodes[k]->evalues);
    }
    
    
    double x,y;
    for(size_t k=0; k<NMODE; k++)
    {
        sprintf(filename,"%s_%i_%i_%i.dat",root,(int)x1,(int)x2,(int)k);
        fptr=fopen(filename,"w");
        for(int n=0; n<=100; n++)
        {
            double theta = atan2(gsl_matrix_get(submodes[k]->evectors,0,1),gsl_matrix_get(submodes[k]->evectors,0,0));
            double Rx = sqrt(gsl_vector_get(submodes[k]->evalues,0));
            double Ry = sqrt(gsl_vector_get(submodes[k]->evalues,1));
            double Cx = gsl_vector_get(submodes[k]->mu,0);
            double Cy = gsl_vector_get(submodes[k]->mu,1);
            
            double angle = n*(2.*M_PI/100.);
            x = 1.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 1.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 2.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 2.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 3.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 3.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }
    
    for(size_t k=0; k<NMODE; k++) free_MVG(submodes[k]);
    free(submodes);
}

void print_model(struct MVG **modes, struct Sample **samples, size_t NMCMC, double logL, double BIC, size_t step)
{
    size_t NP = samples[0]->x->size;
    size_t NMODE = samples[0]->w->size;
    char filename[128];
    
    sprintf(filename,"gmm.dat");
    FILE *fptr = fopen(filename,"a");
    fprintf(fptr,"%i %g %g\n",(int)step, logL, BIC);
    fclose(fptr);
    
    for(size_t m=0; m<NP; m++)
    {
        for(size_t n=m; n<NP; n++)
        {
            /* Get 1D PDFs for plotting */
            if(m==n)
            {
                sprintf(filename,"pdf_%i",(int)step);
                print_1D_pdfs(modes, samples, NMCMC, filename, m);
            }
            
            /* Get Eigenvectors & Eigenvalues of Covariance matrix for plotting */
            else
            {
                sprintf(filename,"contours_%i",(int)step);
                print_2D_contours(modes, NMODE, filename, m, n);
            }
        }
    }
}


void expectation_maximization(struct Sample **samples, struct MVG **modes, size_t NMCMC, double *logL, double *BIC)
{
    size_t NP    = samples[0]->x->size;
    size_t NMODE = samples[0]->p->size;
    
    // aliases to structures
    struct MVG *M = NULL;
    struct Sample *s = NULL;
    
    // helper quantities for building sums etc.
    double norm;
    double mu=0.0;
    double C=0.0;
    double R;
    double p;
    
    /*
     E-step:
     compute probability for each sample to belong to each mode
     */
    
    /* compute p(x|mode) for each sample for each mode */
    for(size_t i=0; i<NMCMC; i++)
    {
        
        norm=0.0;
        s = samples[i];
        
        //loop over modes
        for(size_t k=0; k<NMODE; k++)
        {
            M = modes[k];
            
            //compute p(x|mode)
            p = multivariate_gaussian(s->x, M);
            
            gsl_vector_set(s->p,k,p);
            gsl_vector_set(s->w,k,M->p*p);
            norm += M->p*p;
        }
        gsl_vector_scale(s->w,1./norm);
    }
    
    /* weigh the number of samples in each mode */
    for(size_t k=0; k<NMODE; k++)
    {
        modes[k]->Neff = 0;
        for(size_t i=0; i<NMCMC; i++) modes[k]->Neff += gsl_vector_get(samples[i]->w,k);
        if(modes[k]->Neff < 1.0)
        {
            exit(1);
        }
        modes[k]->p = modes[k]->Neff/(double)NMCMC;
    }
    
    /* check convergence with log likelihood & BIC */
    *logL = log_likelihood(modes, samples, NMCMC, NMODE);
    *BIC = -2.*(*logL) + (double)NMODE*((double)NP*((double)NP+3.)/2. + 1)*log((double)NMCMC);
    printf(" logL = %g,  BIC = %g     ",*logL, *BIC);
    
    
    /*
     M-Step:
     recompute mean, covariance, and relative weight of newly weighted samples
     */
    
    /* compute weighted mean and covariance for each mode */
    for(size_t k=0; k<NMODE; k++)
    {
        M = modes[k];
        
        //get new means
        for(size_t n=0; n<NP; n++)
        {
            mu = 0.0;
            
            //mu is a weighted average for each mode
            for(size_t i=0; i<NMCMC; i++)
                mu += (gsl_vector_get(samples[i]->w,k)/modes[k]->Neff) * gsl_vector_get(samples[i]->x,n);
            
            gsl_vector_set(M->mu,n,mu);
        }
        
        //get new covariance
        for(size_t m=0; m<NP; m++)
        {
            //loop over parameters again (wasteful, not using symmetry
            for(size_t n=0; n<NP; n++)
            {
                C = 0.0;
                
                //loop over samples
                for(size_t i=0; i<NMCMC; i++)
                {
                    double dx_m = gsl_vector_get(samples[i]->x,m) - gsl_vector_get(M->mu,m);
                    double dx_n = gsl_vector_get(samples[i]->x,n) - gsl_vector_get(M->mu,n);
                    
                    C += (gsl_vector_get(samples[i]->w,k)/modes[k]->Neff)*(dx_m)*(dx_n);
                }
                gsl_matrix_set(M->C,m,n,C);
            }
        }
        
        //invert new matrix to evaluate the probabilities
        invert_gsl_matrix(modes[k]->C, modes[k]->Cinv, modes[k]->L, &modes[k]->detC, &R);
    }
}

void GMM_with_EM(struct MVG **modes, struct Sample **samples, size_t NMCMC, size_t NSTEP, gsl_rng *r, double *logL, double *BIC)
{
    size_t NP    = samples[0]->x->size;
    size_t NMODE = samples[0]->p->size;

    /* construct diagonal covariance matrix of full sample variances */
    double x_temp[NMCMC];
    double mean_temp, var_temp;
    for(size_t i=0; i<NP; i++)
    {
        for(size_t n=0; n<NMCMC; n++) x_temp[n] = gsl_vector_get(samples[n]->x,i);
        mean_temp = gsl_stats_mean(x_temp,1,NMCMC);
        var_temp  = gsl_stats_variance_m(x_temp,1,NMCMC, mean_temp);
        
        //set diagonals of C
        for(size_t n=0; n<NMODE; n++) gsl_matrix_set(modes[n]->C,i,i,var_temp);
        
    }
    
    /* place covariance matrices at random draws from the chain file */
    double R; //condition number of matrix
    for(size_t k=0; k<NMODE; k++)
    {
        //pick a sample from the chain to anchor each covariance matrix
        int fair_draw = (int)gsl_ran_flat(r,0,NMCMC);
        for(size_t n=0; n<NP; n++) gsl_vector_set(modes[k]->mu,n,gsl_vector_get(samples[fair_draw]->x,n));
        
        //set priors for each model
        modes[k]->p = (double)1./(double)NMODE;
        
        //get inverset, determinant, etc.
        invert_gsl_matrix(modes[k]->C, modes[k]->Cinv, modes[k]->L, &modes[k]->detC, &R);
    }
    
    /* EM Algorithm for Gaussian Mixture Models */
    size_t step=0;
    double BICmin = 1e60;
    while(step<NSTEP)
    {
        printProgress((double)(step+1)/NSTEP);
        expectation_maximization(samples, modes, NMCMC, logL, BIC);
        if(floor(*BIC) < floor(BICmin))
        {
            BICmin = *BIC;
            step=0;
        }
        step++;
    }
    printf("\n");
}
