//
//  NoiseWrapper.c
//  
//
//  Created by Tyson Littenberg on 2/5/21.
//

#include <string.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#include <LISA.h>

#include <gbmcmc.h>
#include <mbh.h>
#include <Noise.h>

#include "GalacticBinaryWrapper.h"
#include "MBHWrapper.h"
#include "NoiseWrapper.h"

void alloc_noise_data(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, int procID, int nProc)
{
    noise_data->status = 0;
    noise_data->procID = procID;
    noise_data->nProc = nProc;
    noise_data->flags = malloc(sizeof(struct Flags));
    noise_data->chain = malloc(sizeof(struct Chain));
    noise_data->data  = malloc(sizeof(struct Data));
    
    noise_data->orbit = gbmcmc_data->orbit;
    memcpy(noise_data->flags, gbmcmc_data->flags, sizeof(struct Flags));

}



void select_noise_segment(struct Noise *psd_full, struct Data *data, struct Chain *chain, struct Model **model)
{
    double Tobs = data->T;
    int qstart = (int)(psd_full->f[0]*Tobs);
    int q0 = (int)(model[0]->noise[0]->f[0]*Tobs);
    int dq = q0 - qstart;
    
    for(int i=0; i<chain->NC; i++)
    {
        memcpy(model[i]->noise[0]->SnX, psd_full->SnX+dq, data->N*sizeof(double));
        memcpy(model[i]->noise[0]->SnA, psd_full->SnA+dq, data->N*sizeof(double));
        memcpy(model[i]->noise[0]->SnE, psd_full->SnE+dq, data->N*sizeof(double));
    }

}

void setup_noise_data(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, struct MBHData *mbh_data, struct TDI *tdi_full, int procID)
{
    noise_data->data->downsample = gbmcmc_data->data->downsample;
    noise_data->data->Nwave      = gbmcmc_data->data->Nwave;
    
    noise_data->chain->NC      = gbmcmc_data->chain->NC;
    noise_data->data->Nchannel = gbmcmc_data->data->Nchannel;
    noise_data->data->NP       = gbmcmc_data->data->NP;
    strcpy(noise_data->data->format,gbmcmc_data->data->format);

    double T = gbmcmc_data->data->T;
    
    noise_data->data->T = T;
    
    //set noise model to cover gbmcmc segment
    noise_data->data->fmin = gbmcmc_data->data->fmin;
    noise_data->data->fmax = gbmcmc_data->data->fmax;
    
    //set all processes noise models based on max/min ucb segment
    MPI_Bcast(&noise_data->data->fmin, 1, MPI_DOUBLE, gbmcmc_data->procID_min, MPI_COMM_WORLD);
    MPI_Bcast(&noise_data->data->fmax, 1, MPI_DOUBLE, gbmcmc_data->procID_max, MPI_COMM_WORLD);

    //pad noise model
    noise_data->data->fmin /= 1.01;
    noise_data->data->fmax *= 1.01;

    //adjust noise model bandwidth to account for MBHs
    if(mbh_data->NMBH>0)
    {
        //set limits of noise model to cover both models
        noise_data->data->fmin = (mbh_data->data->fmin < noise_data->data->fmin ) ? mbh_data->data->fmin : noise_data->data->fmin;
        noise_data->data->fmax = (mbh_data->data->fmax > noise_data->data->fmax ) ? mbh_data->data->fmax : noise_data->data->fmax;

        //pad noise model even more (MBH bandwidth fluctuates)
        noise_data->data->fmin /= 1.1;
        noise_data->data->fmax *= 1.1;
    }
    //TODO: adjust noise model bandwidth to account for VBs
    
    noise_data->data->N = (int)((noise_data->data->fmax - noise_data->data->fmin)*T);
    
    alloc_data(noise_data->data, noise_data->flags);
    
    noise_data->model = malloc(sizeof(struct SplineModel*)*gbmcmc_data->chain->NC);
    
    //get max and min samples
    noise_data->data->qmin = (int)(noise_data->data->fmin*noise_data->data->T);
    noise_data->data->qmax = noise_data->data->qmin+noise_data->data->N;
    noise_data->data->fmax = (double)noise_data->data->qmax/T;
    
    //store max and min frequency in MBH structure
    mbh_data->data->fmin = noise_data->data->fmin;
    mbh_data->data->fmax = noise_data->data->fmax;
    
    select_frequency_segment(noise_data->data, tdi_full);
    
    /*
     Initialize measured time of model update.
     Used to determine number of steps relative to mbh model
     */
    noise_data->cpu_time = 1.0;

}


void initialize_noise_sampler(struct NoiseData *noise_data)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = noise_data->flags;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize GBMCMC sampler state */
    initialize_noise_state(noise_data);
    
    /* Set sampler counter */
    noise_data->mcmc_step = -flags->NBURN;
    
    /* Store data segment in working directory */
    print_data(data, data->tdi[0], flags, 0);

}

void initialize_noise_state(struct NoiseData *noise_data)
{
    /* Aliases to gbmcmc structures */
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct SplineModel **model = noise_data->model;

    int NC = chain->NC;
    int Nspline = noise_data->nProc+1;
    
    for(int ic=0; ic<NC; ic++)
    {
        model[ic] = malloc(sizeof(struct SplineModel));
        initialize_spline_model(orbit, data, model[ic], Nspline);
    }

    char filename[128];
    sprintf(filename,"%s/initial_spline_points.dat",data->dataDir);
    print_noise_model(model[0]->spline, filename);
    
    sprintf(filename,"%s/interpolated_spline_points.dat",data->dataDir);
    print_noise_model(model[0]->psd, filename);

}

int update_noise_sampler(struct NoiseData *noise_data)
{
    clock_t start = clock();
    
    /* Aliases to gbmcmc structures */
    struct Flags *flags = noise_data->flags;
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct SplineModel **model = noise_data->model;

    int NC = chain->NC;
    
    //For saving the number of threads actually given
    int numThreads;
    
#pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
#pragma omp barrier
        /* The MCMC loop */
                
#pragma omp barrier
        // (parallel) loop over chains
        for(int ic=threadID; ic<NC; ic+=numThreads)
        {
            
            //loop over frequency segments
            struct SplineModel *model_ptr = model[chain->index[ic]];
            
            //update log likelihood (data may have changed)
            model_ptr->logL = noise_log_likelihood(data, model_ptr);

            //evolve fixed dimension sampler
            for(int steps=0; steps<100; steps++)
            {
                noise_spline_model_mcmc(orbit, data, model_ptr, chain, flags, ic);
            }
            
            //evolve trans dimension sampler
            noise_spline_model_rjmcmc(orbit, data, model_ptr, chain, flags, ic);

        }// end (parallel) loop over chains
         
        
        //Next section is single threaded. Every thread must get here before continuing
#pragma omp barrier
        if(threadID==0){
            spline_ptmcmc(model,chain,flags);
            adapt_temperature_ladder(chain, noise_data->mcmc_step+flags->NBURN);
        }
        //Can't continue MCMC until single thread is finished
#pragma omp barrier
        
    }// End of parallelization
    
    print_spline_state(model[chain->index[0]], chain->noiseFile[0], noise_data->mcmc_step);
    
    if(noise_data->mcmc_step>=0 && noise_data->mcmc_step%data->downsample==0 && noise_data->mcmc_step/data->downsample < data->Nwave)
    {
        struct SplineModel *model_ptr = model[chain->index[0]];

        for(int n=0; n<data->N; n++)
        {
            data->S_pow[n][0][0][noise_data->mcmc_step/data->downsample] = model_ptr->psd->SnA[n];
            data->S_pow[n][1][0][noise_data->mcmc_step/data->downsample] = model_ptr->psd->SnE[n];
        }
    }
    
    noise_data->mcmc_step++;

    clock_t stop = clock();
    noise_data->cpu_time = (double)(stop-start);

    return 1;
}
