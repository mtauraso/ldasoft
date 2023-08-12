
/*
 *  Copyright (C) 2019 Kristen Lackeos (MSFC-ST12), Tyson B. Littenberg (MSFC-ST12)
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
 */


/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include <LISA.h>

#include "gbmcmc.h"

// CODE USAGE:
// ./gb_match --match-in1 /Users/klackeos/Desktop/input1.dat --match-in2 /Users/klackeos/Desktop/input2.dat --frac-freq --fmin 0.001249 --samples 512 --duration 62914560

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
    
    FILE *chain_file1;
    FILE *chain_file2;
    int NMAX = 10;   //max number of frequency & time segments
    
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  *data  = malloc(sizeof(struct Data));
    
    
    //   Parse command line and set defaults/flags
    data->t0 = malloc( NMAX * sizeof(double) );
    
    parse(argc,argv,data,orbit,flags,chain,NMAX);
    alloc_data(data, flags);
    data->qmin = (int)(data->fmin*data->T);
    
    /* Load spacecraft ephemerides */
    switch(flags->orbit)
    {
        case 0:
            initialize_analytic_orbit(orbit);
            break;
        case 1:
            initialize_numeric_orbit(orbit);
            break;
        default:
            fprintf(stderr,"unsupported orbit type\n");
            return(1);
            break;
    }
    
    chain_file1 = fopen(flags->matchInfile1,"r");
    chain_file2 = fopen(flags->matchInfile2,"r");
    
    if ( chain_file1 == NULL )
    {
        printf("match-in1 is null\n");
    }
    
    if ( chain_file2 == NULL )
    {
        printf("match-in2 is null\n");
    }
    
    
    //allocate memory for two sources and noise
    struct Source *src1 = NULL;
    src1 = malloc(sizeof(struct Source));
    alloc_source(src1, data->N,2,data->NP);
    
    struct Source *src2 = NULL;
    src2 = malloc(sizeof(struct Source));
    alloc_source(src2, data->N,2,data->NP);
    
    struct Noise *noise = NULL;
    noise = malloc(flags->NT*sizeof(struct Noise));
    alloc_noise(noise, data->N);
    
    
    for(int n=0; n<2*data->N; n++)
    {
        src1->tdi->A[n] = 0.0;
        src1->tdi->E[n] = 0.0;
        src1->tdi->X[n] = 0.0;
        src2->tdi->A[n] = 0.0;
        src2->tdi->E[n] = 0.0;
        src2->tdi->X[n] = 0.0;
    }
    
    //Get noise spectrum for data segment
    for(int n=0; n<data->N; n++)
    {
        double f = data->fmin + (double)(n)/data->T;
        noise->f[n] = f;
        if(strcmp(data->format,"phase")==0)
        {
            noise->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
            noise->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
        }
        else if(strcmp(data->format,"frequency")==0 || strcmp(data->format,"sangria")==0)
        {
            noise->SnA[n] = AEnoise_FF(orbit->L, orbit->fstar, f)/sqrt(2.);
            noise->SnE[n] = AEnoise_FF(orbit->L, orbit->fstar, f)/sqrt(2.);
            noise->SnA[n] += GBnoise_FF(data->T, orbit->fstar, f)/sqrt(2.);
            noise->SnE[n] += GBnoise_FF(data->T, orbit->fstar, f)/sqrt(2.);

        }
        else
        {
            fprintf(stderr,"Unsupported data format %s",data->format);
            exit(1);
        }
    }
    
    
    double max_match;
    double match;
    
    while(!feof(chain_file2))
    {
        for(int n=0; n<2*data->N; n++)
        {
            src2->tdi->A[n] = 0.0;
            src2->tdi->E[n] = 0.0;
        }
        
        scan_source_params(data, src2, chain_file2);
        galactic_binary_alignment(orbit, data, src2);
        galactic_binary(orbit, data->format, data->T, data->t0[0], src2->params, data->NP, src2->tdi->X, src2->tdi->A, src2->tdi->E, src2->BW, 2);

        max_match=-INFINITY;
        while(!feof(chain_file1))
        {

            scan_source_params(data, src1, chain_file1);
            
            if( fabs(src1->params[F0] - src2->params[F0]) < 20.)
            {
                for(int n=0; n<2*data->N; n++)
                {
                    src1->tdi->A[n] = 0.0;
                    src1->tdi->E[n] = 0.0;
                }
                
                //Book-keeping of injection time-frequency volume
                galactic_binary_alignment(orbit, data, src1);
                galactic_binary(orbit, data->format, data->T, data->t0[0], src1->params, data->NP, src1->tdi->X, src1->tdi->A, src1->tdi->E, src1->BW, 2);
                                
                match = waveform_match(src1, src2, noise);
                if(match>max_match) max_match=match;
            }
        }
        printf("%lg %lg %lg\n",max_match,snr(src2,noise),src2->f0);
        rewind(chain_file1);
    }

    fclose(chain_file1);
    fclose(chain_file2);
    
    return 0;
}
