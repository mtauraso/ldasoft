/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

#include <GMM_with_EM.h>

#include <LISA.h>

#include "GalacticBinary.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryCatalog.h"

static inline double galaxy_liklihood(double *x)
{
    double u, rsq, z, s;
    
    z = x[2];
    rsq = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
    u = sqrt(x[0]*x[0]+x[1]*x[1]);
    
    s = 1.0/cosh(z/GALAXY_Zd);
    
    // Note that overall rho0 in density is irrelevant since we are working with ratios of likelihoods in the MCMC
    
    return GALAXY_A*exp(-rsq/(GALAXY_Rb*GALAXY_Rb))+(1.0-GALAXY_A)*exp(-u/GALAXY_Rd)*s*s;
}

static inline double loglike(double *x) {
    return log(galaxy_liklihood(x));
}

// Called by galactic prior mcmc to generate a sample
static inline void generate_galaxy_sample(double *current_sample /*in*/, double *proposed_sample /*out*/, gsl_rng *r) 
{
    double alpha, beta, xx;

    alpha = gsl_rng_uniform(r);
        
    if(alpha > 0.7)  // uniform draw from a galactic bounding box
    {
        _generate_uniform_galaxy_sample(proposed_sample, r);
    }
    else  // Move the current sample a little bit, with a roughly exponential fall-off
    {
        beta = gsl_rng_uniform(r);
        
        if(beta > 0.8)
        {
            xx = 0.1;
        }
        else if (beta > 0.4)
        {
            xx = 0.01;
        }
        else
        {
            xx = 0.001;
        }
        proposed_sample[0] = current_sample[0] + gsl_ran_gaussian(r,xx);
        proposed_sample[1] = current_sample[1] + gsl_ran_gaussian(r,xx);
        proposed_sample[2] = current_sample[2] + 0.1*gsl_ran_gaussian(r,xx);
    }
}

struct Prior * alloc_prior() {
    return calloc(1, sizeof(struct Prior));
}

static inline void calc_sky_buckets(struct Prior *prior, int Nth, int Nph) {
    prior->dcostheta = 2./(double)Nth;
    prior->dphi      = 2.*M_PI/(double)Nph;
    
    prior->ncostheta = Nth;
    prior->nphi      = Nph;
}

static inline void alloc_sky_prior(struct Prior * prior, int Nth, int Nph) {
    // Note that calloc's practice of setting memory to 0
    // also sets IEEE-754 double precision floating point values to +0.0
    prior->skyhist = (double*)calloc((Nth*Nph),sizeof(double));
    calc_sky_buckets(prior, Nth, Nph);
    prior->skymaxp = 0.0;
}

static inline void alloc_volume_prior(struct Prior *prior, struct Flags *flags, int Nth, int Nph, int Nr) {
    prior->spherehist = (double*)calloc((Nth*Nph*Nr), sizeof(double));
    prior->spheremaxp = 0.0;

    calc_sky_buckets(prior, Nth, Nph);
    prior->nr = Nr;
    prior->dr = (GALAXY_BS_R)/Nr;

    prior->fdotastroPrior = flags->fdotastroPrior;
}

// Takes in a sky direciton (costheta, phi) returns the corresponding indicies for costheta and phi 
// suitable for either the sky prior or the volume prior.
static inline void sky_direction_to_sky_bucket(struct Prior *prior, double costheta, double phi, int *ith, int *iph) {
    int Nth = prior->ncostheta;
    int Nph = prior->nphi;
    *ith = (int)(0.5*(1.0+costheta)*(double)(Nth));
    //*iph = (int)((2*M_PI-phi)/(2.0*M_PI)*(double)(Nph));
    *iph = (int)((phi)/(2.0*M_PI)*(double)(Nph));
    
    // Error checking
    if(*ith < 0 || *ith > Nth -1) printf("Out of bounds in Theta direction: %d %d\n", *ith, *iph);
    if(*iph < 0 || *iph > Nph -1) printf("Out of bounds in Phi direction: %d %d\n", *ith, *iph);
}

// Takes sky-distance representation (as it is in params) and returns
// the single index in the spherical volume bucket corresponding to that point
static inline int sky_distance_to_sphere_index(struct Prior *prior, double costheta, double phi, double r_ec) {
    int Nph = prior->nphi;
    int Nr = prior->nr;

    int ith, iph, ir;
    sky_direction_to_sky_bucket(prior, costheta, phi, &ith, &iph);
    ir = (int)(r_ec/prior->dr);

    if(ir < 0 || ir > Nr -1) {
        printf("Out of bounds in R direction: %d %d %d\n", ith, iph, ir);
    }

    return ith*(Nph*Nr) + iph*(Nr) + ir;
}

// Called by galactic prior generator mcmc to sample the chain to generate the sky prior
static inline void sample_sky_prior(struct Prior *prior, double *chain_sample)
{
    int ith, iph;
    double phi, theta, r_ec;
    int Nph = prior->nphi;

    galactocentric_to_sky_distance(chain_sample, &phi, &theta, &r_ec);

    // sin(theta) here because the "theta" returned by galactocentric_to_sky_distance is 
    // an elevation angle measured from the celestial equator with range [-pi/2, pi/2]
    // This function expects a coordinate which is cos("theta" + pi/2) with the range [-1,1]
    // as is typical of params, e.g. The cosine of an angle measured positively from the north pole
    sky_direction_to_sky_bucket(prior, sin(theta), phi, &ith, &iph);

    prior->skyhist[ith*Nph+iph] += 1.0;
}

// Called by the galactic prior generator to sample the chain and generate the 3d volumetric galaxy prior
static inline void sample_volume_prior(struct Prior *prior, double *chain_sample) {
    double phi, theta, r_ec;
    galactocentric_to_sky_distance(chain_sample, &phi, &theta, &r_ec);
    
    // sin(theta) here because the "theta" returned by galactocentric_to_sky_distance is 
    // an elevation angle measured from the celestial equator with range [-pi/2, pi/2]
    // This function expects a coordinate which is cos("theta" + pi/2) with the range [-1,1]
    // as is typical of params, e.g. The cosine of an angle measured positively from the north pole
    int i = sky_distance_to_sphere_index(prior, sin(theta), phi, r_ec);

    // Index, add one to the count
    prior->spherehist[i] += 1.0;
}

static inline void sphere_index_to_coords(struct Prior *prior, int index, double *costheta, double *phi, double *r_ec ) {

    // Convert our single index into indicies along costheta, phi, r directions of our 3d array.
    int costheta_idx = (index / (prior->nphi*prior->nr));
    int phi_idx =      (index % (prior->nphi*prior->nr))/ prior->nr;
    int r_idx =        (index % (prior->nphi*prior->nr))% prior->nr;

    // add 0.5 to each idx before converting to get the center of each dV region
    *costheta = ((costheta_idx+0.5)*prior->dcostheta) - 1.0;
    *phi =      ((phi_idx+0.5)*prior->dphi);
    *r_ec =     ((r_idx+0.5)*prior->dr);
}

// Called as the last step in generation of the sky prior.
// Histogram values are converted from counts to log(probability / solid angle)
// 10% of total probability is subtracted from counts and distributed evenly across all bins (set with uni)
// Sets maximum histogram value for sky prior
void convert_sky_prior(struct Prior *prior, int cnt) {
    double dOmega, xx, yy, zz;
    int iph, ith;
    int Nth = prior->ncostheta;
    int Nph = prior->nphi;

    // Solid angle subtended by an angular bin
    dOmega = 4.0*M_PI/(double)(Nth*Nph);
    
    double uni = 0.1;
    //fprintf(stderr,"\n   HACK:  setup_galaxy_prior() uni=%g\n",uni);

    yy = (1.0-uni)/(double)(cnt);
    zz = uni/(double)(Nth*Nph);
    
    yy /= dOmega;
    zz /= dOmega;
    
    for(ith=0; ith< Nth; ith++)
    {
        for(iph=0; iph< Nph; iph++)
        {
            xx = yy*prior->skyhist[ith*Nph+iph];
            prior->skyhist[ith*Nph+iph] = log(xx + zz);
            if(prior->skyhist[ith*Nph+iph]>prior->skymaxp) prior->skymaxp = prior->skyhist[ith*Nph+iph];
        }
    }
}

// xcxc: Ought this be something more like a plummer potential
//       to reduce support right on top of earth?
//
// We allocate additional probability to the volume prior in order to provide
// support for sky locations outside the galaxy. Volume prior uses the below 
// profile function to allocate this probability. 
//
// This function need not be normalized, and is parameterized in terms of r_ec
// so it will be anisotropic over the sky (with the exception of any effects 
// introduced by volume prior bounding box alignment)
// 
// At long distance this function should fall off at least as rapidly as 1/r^2 
// to prevent the prior from biasing toward more distant sources.
static inline double volume_prior_uniform_contribution(double r_ec){
    // something more plummer-like
    double a = 10; //kpc
    return 1/(r_ec*r_ec + a*a);

    // 1/r^2
    //return 1/(r_ec*r_ec);

    // log normal profile

}

// Calculate the triple integral of our uniform contribution function
// added up across every box in the volume.
double volume_prior_uniform_normalization(struct Prior *prior) {
    
    int num_buckets = prior->ncostheta*prior->nphi*prior->nr;
    double integral = 0.0;

    fprintf(stdout,"Normalizing Volume Prior Uniform Contribution\n");
    for(int i=0; i< num_buckets; i++) {
        double phi, costheta, r_ec;
        sphere_index_to_coords(prior, i, &costheta, &phi, &r_ec);

        // volume_prior_uniform_contribution has units probability * kpc^-3
        // Therefore we must sum it over the spherical volume element
        // to get a total probability.
        double dVol = r_ec * r_ec * prior->dr * prior->dcostheta * prior->dphi;
        integral += volume_prior_uniform_contribution(r_ec) * dVol;
    }

    return 1.0/integral;
}

// Return 1/(integral) where integral is exp(loglike(x))
//  for all X in our spherical volume
double galaxy_liklihood_normalization(struct Prior *prior) {
    int num_buckets = prior->ncostheta*prior->nphi*prior->nr;
    double integral = 0.0;

    fprintf(stdout,"Normalizing Volume Prior Galaxy Contribution\n");
    for(int i=0; i< num_buckets; i++) {
        double x[3];
        double phi, costheta, r_ec;
        sphere_index_to_coords(prior, i, &costheta, &phi, &r_ec);
        sky_distance_to_galactocentric(x, phi, M_PI/2.0 - acos(costheta), r_ec);

        // Units of galaxy_liklihood are probability * kpc^-3 
        // So we must use the spherical volume element when summing them.
        double dVol = r_ec * r_ec * prior->dr * prior->dcostheta * prior->dphi;

        integral += galaxy_liklihood(x) * dVol;
    }
    return 1.0/integral;
}

void convert_volume_prior(struct Prior *prior, int cnt) {
    int num_buckets = prior->ncostheta*prior->nphi*prior->nr;

    // Amount of total probabilty to redistribute uniformly across volume
    double uni = 0.1;

    // Normalizing factor for converting count-> probability but with unitary 
    // contribution subtracted out.
    double count_normalization = (1.0 - uni)/(double)(cnt);

    // normalization for the uniform-over-sky-angle contribution
    double uniform_normalization = uni * volume_prior_uniform_normalization(prior);

    for(int i = 0; i < num_buckets; i++) {
        double phi, costheta, r_ec;
        sphere_index_to_coords(prior, i, &costheta, &phi, &r_ec);
        double dVol = r_ec * r_ec * prior->dr * prior->dcostheta * prior->dphi;
        double uniform_contribution = uniform_normalization*volume_prior_uniform_contribution(r_ec)/dVol;
        double mcmc_contribution = count_normalization*prior->spherehist[i]/dVol;

        prior->spherehist[i] = log(uniform_contribution + mcmc_contribution);
        if(prior->spherehist[i] > prior->spheremaxp) prior->spheremaxp = prior->spherehist[i];
    }
}


void compute_volume_prior(struct Prior *prior) {
    int num_buckets = prior->ncostheta*prior->nphi*prior->nr;

    // Amount of total probabilty to redistribute uniformly across volume
    double uni = 0.1;

    // normalization for the uniform-over-sky-angle contribution
    double uniform_normalization = uni * volume_prior_uniform_normalization(prior);
    double galaxy_normalization = (1.0 - uni) * galaxy_liklihood_normalization(prior);

    // todo: Change memory layout so we iterate phi/theta last
    //       with r as the major coordinate, then make this a double loop
    //       outer loop iterates r, calculates dVol, inner loop does the rest.
    //       Should speed this up quite a bit from reduced calculations and 
    //       memory locality
    fprintf(stdout,"Computing spherical galaxy prior\n");
    for(int i = 0; i < num_buckets; i++) {
        if(i%(num_buckets/100)==0) printProgress((double)i/(double)(num_buckets));

        double phi, costheta, r_ec;
        double x[3];
        sphere_index_to_coords(prior, i, &costheta, &phi, &r_ec);

        // This function expects an angle measured as elevation from ecliptic equator
        sky_distance_to_galactocentric(x, phi, M_PI/2.0 - acos(costheta), r_ec);

        double uniform_contribution = uniform_normalization*volume_prior_uniform_contribution(r_ec);
        double galaxy_contribution = galaxy_normalization*galaxy_liklihood(x);

        // This has units of log(probability * kpc^-3)
        prior->spherehist[i] = log(uniform_contribution + galaxy_contribution);
        // for debugging just one part of the probability distribution
        //prior->spherehist[i] = log(uniform_contribution);
        //prior->spherehist[i] = log(galaxy_contribution);
        
        if(prior->spherehist[i] > prior->spheremaxp) prior->spheremaxp = prior->spherehist[i];
    }
    printProgress(1.0);
    fprintf(stdout,"\n");
}

// Called in verbose mode to dump the sky prior to a file
// Output sky prior in verbsoe mode
void dump_sky_prior(struct Prior *prior)
{
    int ith, iph;
    int Nth = prior->ncostheta;
    int Nph = prior->nphi;
    double xx, yy;
    FILE *fptr = fopen("skyprior.dat", "w");
    
    for(ith=0; ith< Nth; ith++)
    {
        xx = -1.0+2.0*((double)(ith)+0.5)/(double)(Nth);
        for(iph=0; iph< Nph; iph++)
        {
            yy = 2.0*M_PI*((double)(iph)+0.5)/(double)(Nph);
            fprintf(fptr,"%e %e %e\n", yy, xx, prior->skyhist[ith*Nph+iph]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

void dump_volume_prior(struct Prior *prior)
{
    int num_buckets = prior->ncostheta*prior->nphi*prior->nr;
    FILE *fptr = fopen("3dsphereprior.dat", "w");
    fprintf(fptr, "# costheta phi r_ec log_density \n");

    FILE *fptrcart = fopen("3dvolumeprior.dat", "w");
    fprintf(fptrcart, "# x y z log_density \n");

    fprintf(stdout,"Writing spherical galaxy prior to disk.\n");

    for(int i = 0; i < num_buckets; i++) {
        if(i%(num_buckets/100)==0) printProgress((double)i/(double)(num_buckets));
        double phi, costheta, r_ec;
        double x[3];
        sphere_index_to_coords(prior, i, &costheta, &phi, &r_ec);
        // This function expects an angle measured as elevation from ecliptic
        sky_distance_to_galactocentric(x, phi, M_PI/2.0 - acos(costheta), r_ec);
        fprintf(fptr, "%e %e %e %e\n", costheta, phi, r_ec, prior->spherehist[i]);
        fprintf(fptrcart, "%e %e %e %e\n", x[0], x[1], x[2], prior->spherehist[i]);
    }
    printProgress(1.0);
    fclose(fptr);
    fclose(fptrcart);
}

// assumes flags->volumePrior = true;
void set_volume_prior(struct Flags *flags, struct Prior *prior) 
{
    int Nth = 200;  // bins in cos theta
    int Nph = 200;  // bins in phi
    int Nr = 200;  // bins in r_ec (if applicable)

    if(!flags->quiet)
    {
        fprintf(stdout,"\n============ Galaxy model spherical volumetric prior ============\n");
        fprintf(stdout,"Analytic galaxy model\n");
        fprintf(stdout,"   Distance to GC  = %g kpc\n",GALAXY_RGC);
        fprintf(stdout,"   Disk Radius     = %g kpc\n",GALAXY_Rd);
        fprintf(stdout,"   Disk Height     = %g kpc\n",GALAXY_Zd);
        fprintf(stdout,"   Bulge Radius    = %g kpc\n",GALAXY_Rb);
        fprintf(stdout,"   Bulge Fraction  = %g\n",    GALAXY_A);
        fprintf(stdout,"   Radus around earth (kpc) = %g\n", GALAXY_BS_R);
        fprintf(stdout,"   Bin counts (costheta, phi, r) = %i, %i, %i\n", Nth, Nph, Nr);
    }

    alloc_volume_prior(prior, flags, Nth, Nph, Nr);

    compute_volume_prior(prior);

    if(flags->verbose)
    {
        dump_volume_prior(prior);
    }
}

void set_galaxy_prior(struct Flags *flags, struct Prior *prior)
{
    if(!flags->quiet)
    {
        if(flags->galaxyPrior) fprintf(stdout,"\n============ Galaxy model sky prior ============\n");
        if(flags->volumePrior) fprintf(stdout,"\n========= Galaxy model volumetric prior ========\n");

        fprintf(stdout,"Monte carlo over galaxy model\n");
        fprintf(stdout,"   Distance to GC  = %g kpc\n",GALAXY_RGC);
        fprintf(stdout,"   Disk Radius     = %g kpc\n",GALAXY_Rd);
        fprintf(stdout,"   Disk Height     = %g kpc\n",GALAXY_Zd);
        fprintf(stdout,"   Bulge Radius    = %g kpc\n",GALAXY_Rb);
        fprintf(stdout,"   Bulge Fraction  = %g\n",    GALAXY_A);
        fprintf(stdout,"   Bounding Volume = %g x %g x %g kpc\n", GALAXY_BB_X, GALAXY_BB_Y, GALAXY_BB_Z);
    }
    double *x, *y;  // current and proposed parameters
    int D = 3;  // number of parameters
    int Nth = 200;  // bins in cos theta
    int Nph = 200;  // bins in phi
    int Nr = 200;  // bins in r_ec (if applicable)
    int MCMC=100000000;
    int j;
    int cnt;
    double H;
    double logLx, logLy;
    double beta;
    int mc;
    //  FILE *chain = NULL;
    
    if(flags->debug)
    {
        Nth /= 10;
        Nph /= 10;
        MCMC/=10;
    }
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    x =  (double*)calloc(D,sizeof(double));
    y =  (double*)calloc(D,sizeof(double));
    
    if(flags->galaxyPrior) alloc_sky_prior(prior, Nth, Nph);
    if(flags->volumePrior) alloc_volume_prior(prior, flags, Nth, Nph, Nr);
    
    // starting values for parameters
    x[0] = 0.5;
    x[1] = 0.4;
    x[2] = 0.1;
    
    logLx = loglike(x);
    
    cnt = 0;
    
    for(mc=0; mc<MCMC; mc++)
    {
        if(mc%(MCMC/100)==0) printProgress((double)mc/(double)MCMC);

        generate_galaxy_sample(x, y, r);

        logLy = loglike(y);
        
        H = logLy - logLx;

        beta = gsl_rng_uniform(r);
        beta = log(beta);
        
        if(H > beta)
        {
            for(j=0; j< D; j++) x[j] = y[j];
            logLx = logLy;
        }
        
        if(mc%100 == 0)
        {
            // TODO flags -> --3d-galaxy-prior, only do one of these sampling processes.
            if(flags->galaxyPrior) sample_sky_prior(prior, x);
            if(flags->volumePrior) sample_volume_prior(prior, x);
            cnt++;
        }
    }

    if(flags->galaxyPrior) convert_sky_prior(prior, cnt);
    if(flags->volumePrior) convert_volume_prior(prior, cnt);

    // Output sky prior in verbsoe mode
    if(flags->verbose)
    {
        if(flags->galaxyPrior) dump_sky_prior(prior);
        if(flags->volumePrior) dump_volume_prior(prior);
    }
    
    free(x);
    free(y);
    gsl_rng_free (r);
    if(!flags->quiet)fprintf(stdout,"\n================================================\n\n");
    fflush(stdout);
    
}

// This arose from a curve examination of the fdot_astro of
// interacting LDC injected sources 
// Some educated guesses were made about M_c
//
// The calculations showing 1/(x^2 + a^2) is a reasonable profile
// across a few choices of M_c are shown in scripts/fdot_prior.ipynb
//
# define DFDTASTRO_A (8e-21) // Hz^2

// singleton to store the norm.
// There's not a great place to put this given its dependent on 
// the limits on dfdt stored in model->prior, but thematically is best
// put into struct Prior.
static double DFDT_ASTRO_PRIOR_NORM = NAN;

static inline double dfdtastro_prior_norm(struct Data * data, struct Model *model) {
    assert(is_param(DFDTASTRO));

    // if its not Nan, we already ran once and set it.
    if(DFDT_ASTRO_PRIOR_NORM == DFDT_ASTRO_PRIOR_NORM) {
        return DFDT_ASTRO_PRIOR_NORM;
    }

    // convert a to data dependent frequency units
    double a = DFDTASTRO_A*data->T*data->T;

    double fdotastromin = model->prior[DFDTASTRO][0];
    double fdotastromax = model->prior[DFDTASTRO][1];

    // This is the log of total volume under the curve
    // (x^2 + a^2)^-1 between fdotastromin and fdotastromax
    // in the same units as params.
    //
    // log{ a^-1 * (atan(fdotmax/a) - atan(fdotastromin/a)) }
    DFDT_ASTRO_PRIOR_NORM = -log(a) + log(atan(fdotastromax/a) - atan(fdotastromin/a));

    return DFDT_ASTRO_PRIOR_NORM;
}

double evaluate_dfdtastro_prior(struct Prior *prior, struct Data * data, struct Model* model, double *params) {
    assert(is_param(DFDTASTRO));

    if(prior->fdotastroPrior) {
        fprintf(stderr, "Currently fdotastro nonuniform prior has not been checked for units, and should not be used.\n");
        exit(1);

        // TODO: Need to think through given arguments about units of top/bottom of acceptance
        //       ratio, whether the non-uniform case (prior->fdtoastroPrior = true) works correctly
        //       here.
        
        double a = DFDTASTRO_A*data->T*data->T; // convert a to params units

        double a_sq = a*a;
        double dfdt_sq = params[DFDTASTRO]*params[DFDTASTRO];

        // log{ (total_area)^-1 * (dfdt^2 + a^2)^-1 }
        return - dfdtastro_prior_norm(data, model) - log(dfdt_sq + a_sq);
    } else {
        return -model->logPriorVolume[DFDTASTRO];
    }
}

void set_uniform_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose)
{
    /*
     params[F0] = source->f0*T;
     params[COSTHETA] = source->costheta;
     params[PHI] = source->phi;
     params[AMP] = log(source->amp);
     params[COSI] = source->cosi;
     params[PSI] = source->psi;
     params[PHI0] = source->phi0;
     params[DFDT] = source->dfdt*T*T;
     
     See map_array_to_params and map_params_to_array for other parameters.
     */
    
    
    //TODO:  make t0 a parameter
    for(int i=0; i<model->NT; i++)
    {
        model->t0[i] = data->t0[i];
        model->t0_min[i] = data->t0[i] - 20;
        model->t0_max[i] = data->t0[i] + 20;
    }
    
    //TODO: assign priors by parameter name, use mapper to get into vector (more robust to changes)
    
    //frequency bin
    model->prior[F0][0] = data->qmin;
    model->prior[F0][1] = data->qmax;
    
    //colatitude
    model->prior[COSTHETA][0] = -1.0;
    model->prior[COSTHETA][1] =  1.0;
    
    //longitude
    model->prior[PHI][0] = 0.0;
    model->prior[PHI][1] = PI2;
    
    if(is_param(AMP)) {
        //log amplitude
        model->prior[AMP][0] = LOG_AMP_MIN;
        model->prior[AMP][1] = LOG_AMP_MAX;
    }
    
    //cos inclination
    model->prior[COSI][0] = -1.0;
    model->prior[COSI][1] =  1.0;
    
    //polarization
    model->prior[PSI][0] = 0.0;
    model->prior[PSI][1] = M_PI;
    
    //phase
    model->prior[PHI0][0] = 0.0;
    model->prior[PHI0][1] = PI2;
    
    double Mcmin = 0.15;
    double Mcmax = 1.00;
    if(is_param(MC)) {
        model->prior[MC][0] = Mcmin;
        model->prior[MC][1] = Mcmax;
    }

    if(is_param(DIST)) {
        model->prior[DIST][0] = 0;
        model->prior[DIST][1] = GALAXY_BS_R;
    }

    //fdot (bins/Tobs)
    
    /* frequency derivative priors are a little trickier...*/
    double fmin = model->prior[F0][0]/data->T;
    double fmax = model->prior[F0][1]/data->T;
    
    /* emprical envelope functions from Gijs' MLDC catalog */
    double fdotmin = -0.000005*pow(fmin,(13./3.));
    double fdotmax = 0.0000008*pow(fmax,(11./3.));
    
    /* unphysically broad priors
    double fdotmin = -pow(fmin,(13./3.));
    double fdotmax = pow(fmax,(13./3.)); */

    /* use prior on chirp mass to convert to priors on frequency evolution */
    /* driven by gravitational waves alone*/
    double fdotgrmin = galactic_binary_fdot(Mcmin, fmin);
    double fdotgrmax = galactic_binary_fdot(Mcmax, fmax); 
    
    // Ensure our prior on fdotastro allows any mc parameter to hit the full 
    // range of [fdotmin, fdotmax]. This the least constraining version of the
    // fdot astro prior
    // when calculating fdot = fdot_gr + fdot_astro
    //double fdotastromin = fdotmin - fdotgrmax;
    //double fdotastromax = fdotmax - fdotgrmin;

    //Ensure our prior on fdotastro has support only in a narrow range around zero
    // that still allows fdotmin and fdotmax to be hit at parameter extremes of mc
    double fdotastromin = fdotmin - fdotgrmin;
    double fdotastromax = fdotmax - fdotgrmax;
    
    
    if(flags->detached)
    {
        fdotmin = fdotgrmin;
        fdotmax = fdotgrmax;
    }
    
    double fddotmin = 11.0/3.0*fdotmin*fdotmin/fmax;
    double fddotmax = 11.0/3.0*fdotmax*fdotmax/fmin;
    
    if(!flags->detached)
    {
        fddotmin = -fddotmax;
    }


    if(verbose && !flags->quiet)
    {
        fprintf(stdout,"\n============== PRIORS ==============\n");
        if(flags->detached)fprintf(stdout,"  Assuming detached binary, Mchirp = [0.15,1]\n");
        
        if(is_param(DFDT))
        fprintf(stdout,"  p(fdot)      = U[%g,%g]\n",fdotmin,fdotmax);

        if(is_param(DFDTASTRO) && flags->fdotastroPrior)
        fprintf(stdout,"  p(fdotastro) = Plummer[%g,%g,%g]\n",fdotastromin,fdotastromax,DFDTASTRO_A*data->T*data->T);
        else if(is_param(DFDTASTRO))
        fprintf(stdout,"  p(fdotastro) = U[%g,%g]\n",fdotastromin,fdotastromax);
        
        if(is_param(D2FDT2))
        fprintf(stdout,"  p(fddot)     = U[%g,%g]\n",fddotmin,fddotmax);

        if(is_param(AMP)) 
        fprintf(stdout,"  p(lnA)        = U[%g,%g]\n",model->prior[AMP][0],model->prior[AMP][1]);

        if(is_param(MC)) 
        fprintf(stdout,"  p(Mc)         = U[%g,%g]\n",model->prior[MC][0],model->prior[MC][1]);
        fprintf(stdout,"====================================\n\n");
    }
    
    if(is_param(DFDT))
    {
        model->prior[DFDT][0] = fdotmin*data->T*data->T;
        model->prior[DFDT][1] = fdotmax*data->T*data->T;
    }

    if(is_param(DFDTASTRO)) {
        if(flags->detached) {
            // in detached mode fdot from astrophysical sources is zero
            model->prior[DFDTASTRO][0] = 0.0;
            model->prior[DFDTASTRO][1] = 0.0;
        } else {        
            model->prior[DFDTASTRO][0] = fdotastromin*data->T*data->T;
            model->prior[DFDTASTRO][1] = fdotastromax*data->T*data->T;
        }
    }

    if(is_param(D2FDT2))
    {
        model->prior[D2FDT2][0] = fddotmin*data->T*data->T*data->T;
        model->prior[D2FDT2][1] = fddotmax*data->T*data->T*data->T;
    }
    
    /*
     
     Learn to parse prior files with flexible format and
     reset uniform priors accordingly (doing the naive thing
     of setting a uniform distribution over the ~90% credible
     interval for example.  Will need to expand to at least
     gaussian priors.
     
     What to do about intervals so small that they are
     effectively delta functions?  Might anger some proposals...
     
     GAIA distance accuracy < 20%
     
     */
    if(flags->emPrior)
    {
        FILE *priorFile = fopen(flags->pdfFile,"r");
        char name[16];
        double min,max;
        bool parse_error = false;
        
        while(fscanf(priorFile,"%s %lg %lg",name,&min,&max) != EOF)
        {
            
            if(strcmp("f0",name) == 0)
            {
                model->prior[F0][0] = min*data->T;
                model->prior[F0][1] = max*data->T;
            }
            
            else if(strcmp("costheta",name) == 0)
            {
                model->prior[COSTHETA][0] = min;
                model->prior[COSTHETA][1] = max;
            }
            
            else if(strcmp("phi",name) == 0)
            {
                model->prior[PHI][0] = min;
                model->prior[PHI][1] = max;
            }

            else if(strcmp("amp",name) == 0 && is_param(AMP))
            {
                model->prior[AMP][0] = log(min);
                model->prior[AMP][1] = log(max);
            }
            else if(strcmp("amp",name) == 0)
            {
                fprintf(stderr,"amp specified in EM prior file, but in the given command line that is not a search parameter\n");
                parse_error = true;
            }
            
            else if(strcmp("cosi",name) == 0)
            {
                model->prior[COSI][0] = min;
                model->prior[COSI][1] = max;
            }
            
            else if(strcmp("psi",name) == 0)
            {
                model->prior[PSI][0] = min;
                model->prior[PSI][1] = max;
            }
            
            else if(strcmp("phi0",name) == 0)
            {
                model->prior[PHI0][0] = min;
                model->prior[PHI0][1] = max;
            }

            else if(strcmp("dfdt",name) == 0 && is_param(DFDT))
            {
                model->prior[DFDT][0] = min*data->T*data->T;
                model->prior[DFDT][1] = max*data->T*data->T;
            }
            else if(strcmp("dfdt",name) == 0) 
            {
                fprintf(stderr,"dfdt specified in EM prior file, but in the given command line that is not a search parameter\n");
                parse_error = true;
            }

            else if(strcmp("d2fdt2",name) == 0 && is_param(D2FDT2))
            {
                model->prior[D2FDT2][0] = min*data->T*data->T*data->T;
                model->prior[D2FDT2][1] = max*data->T*data->T*data->T;
            }
            else if(strcmp("d2fdt2",name) == 0)
            {
                fprintf(stderr,"d2fdt2 specified in EM prior file, but in the given command line that is not a search parameter\n");
                parse_error = true;
            }

            else if(strcmp("dist",name) == 0 && is_param(DIST))
            {
                model->prior[DIST][0] = min;
                model->prior[DIST][1] = max;
            }
            else if(strcmp("dist",name) == 0)
            {
                fprintf(stderr,"dist specified in EM prior file, but in the given command line that is not a search parameter\n");
                parse_error = true;
            }

            else if(strcmp("mchirp",name) == 0 && is_param(MC))
            {
                model->prior[DIST][0] = min;
                model->prior[DIST][1] = max;
            }
            else if(strcmp("mchirp",name) == 0)
            {
                fprintf(stderr,"mchirp specified in EM prior file, but in the given command line that is not a search parameter\n");
                parse_error = true;
            }

            else if(strcmp("dfdtastro",name) == 0 && is_param(DFDTASTRO))
            {
                model->prior[DFDTASTRO][0] = min*data->T*data->T;
                model->prior[DFDTASTRO][1] = max*data->T*data->T;
            }
            else if(strcmp("dfdtastro",name) == 0) 
            {
                // TODO: xcxc
                // If mchirp, dfdt, f0 have given priors we can technically calculate a range as is done
                // for the in-built prior rather than erroring out.
                fprintf(stderr,"dfdtastro specified in EM prior file, but in the given command line that is not a search parameter\n");
                parse_error = true;
            }
            
            else
            {
                fprintf(stdout, "unrecognized parameter in prior file: %s\n",name);
                parse_error = true;
            }      
        }

        fclose(priorFile);
        if(parse_error) exit(1);
    }
    
    //set prior volume
    for(int n=0; n<data->NP; n++) {
        // Log Prior volume for distance is not constant over distance and is handled in distance_draw_logP
        if(n == DIST) model->logPriorVolume[n] = NAN;
        else model->logPriorVolume[n] = log(model->prior[n][1]-model->prior[n][0]);
    }
    
}

enum SkyPriorMode get_sky_prior_mode(struct Flags *flags) {
    if(flags->galaxyPrior)      return galaxyPrior;
    else if(flags->volumePrior) return volumePrior;                
    return uniformPrior;
}

int check_range(double *params, double **uniform_prior, int NP)
{
    //nan check
    for(int n=0; n<NP; n++) if(params[n]!=params[n]) return 1;
    
    //frequency bin (uniform)
    if(params[F0]<uniform_prior[F0][0] || params[F0]>uniform_prior[F0][1]) return 1;
    
    //cosine co-latitude
    if(params[COSTHETA]<uniform_prior[COSTHETA][0] || params[COSTHETA]>uniform_prior[COSTHETA][1]) return 1;

    //longitude
    if(params[PHI]<uniform_prior[PHI][0] || params[PHI]>=uniform_prior[PHI][1])
    {
        params[PHI] = atan2(sin(params[PHI]),cos(params[PHI]));
        if(params[PHI] < 0.0) params[PHI] += PI2;
    }

    //cosine inclination
    if(params[COSI]<uniform_prior[COSI][0] || params[COSI]>uniform_prior[COSI][1]) return 1;
    
    //polarization
    if(params[PSI]<uniform_prior[PSI][0] || params[PSI]>uniform_prior[PSI][1]) params[PSI] = fmod(params[PSI], M_PI);
    // fmod of a negative first arg returns a negative, this should fix with a small number of iterations
    while(params[PSI]<uniform_prior[PSI][0]) params[PSI] += M_PI;
    if(params[PSI]<uniform_prior[PSI][0] || params[PSI]>uniform_prior[PSI][1]) return 1;

    //phase
    if(params[PHI0]<uniform_prior[PHI0][0] || params[PHI0]>uniform_prior[PHI0][1]) params[PHI0] = fmod(params[PHI0], PI2);
    // fmod of a negative first arg returns a negative, this should fix value with a small number of iterations
    while(params[PHI0]<uniform_prior[PHI0][0]) params[PHI0] += PI2;
    if(params[PHI0]<uniform_prior[PHI0][0] || params[PHI0]>uniform_prior[PHI0][1]) return 1;
    
    //fdot (bins/Tobs)
    if(is_param(DFDT) && (params[DFDT]<uniform_prior[DFDT][0] || params[DFDT]>uniform_prior[DFDT][1])) return 1;
    if(is_param(DFDTASTRO) && (params[DFDTASTRO]<uniform_prior[DFDTASTRO][0] || params[DFDTASTRO]>uniform_prior[DFDTASTRO][1])) return 1;
    
    //fddot
    if(is_param(D2FDT2) && (params[D2FDT2]<uniform_prior[D2FDT2][0] || params[D2FDT2]>uniform_prior[D2FDT2][1])) return 1;

    //mchirp
    if(is_param(MC) && (params[MC]<uniform_prior[MC][0] || params[MC]>uniform_prior[MC][1])) return 1;

    return 0;
}

void set_gmm_prior(struct Flags *flags, struct Data *data, struct Prior *prior)
{
    //get size of full catalog
    int N = data->catalog->N;
    
    //allocate gmm to include the full catalog
    prior->gmm = malloc(sizeof(struct GMM));
    prior->gmm->NP = data->catalog->entry[0]->gmm->NP;
    prior->gmm->NMODE = 0;
    for(size_t n=0; n<N; n++) prior->gmm->NMODE += data->catalog->entry[n]->gmm->NMODE;
    prior->gmm->modes = malloc(prior->gmm->NMODE*sizeof(struct MVG *));
    
    for(size_t n=0; n<prior->gmm->NMODE; n++)
    {
        prior->gmm->modes[n] = malloc(sizeof(struct MVG));
        
        alloc_MVG(prior->gmm->modes[n], (size_t)prior->gmm->NP);
    }

    //combine modes into one GMM
    size_t m=0;
    for(size_t n=0; n<N; n++)
    {
        for(size_t i=0; i<data->catalog->entry[n]->gmm->NMODE; i++)
        {
            copy_MVG(data->catalog->entry[n]->gmm->modes[i],prior->gmm->modes[m]);

            //(clumsily) renormalize modes
            prior->gmm->modes[m]->p /= (double)N;
            
            m++;
        }
    }
    
    //prior->gmm = data->catalog->entry[0]->gmm;
}

void free_prior(struct Prior *prior) {
    if(prior->skyhist) free(prior->skyhist);
    if(prior->spherehist) free(prior->spherehist);
    if(prior->gmm) {
        for(size_t n=0; n<prior->gmm->NMODE; n++) free_MVG(prior->gmm->modes[n]);
        free(prior->gmm->modes);
        free(prior->gmm);
    }
    free(prior);
}

double evaluate_gmm_prior(struct Data *data, struct GMM *gmm, double *params)
{
    size_t NP = data->NP;
    gsl_vector *x = gsl_vector_alloc(NP);
    
    /* pointers to GMM contents */
    struct MVG **modes = gmm->modes;
    size_t NMODES = gmm->NMODE;
    
    //pack parameters into gsl_vector with correct units
    struct Source *source = malloc(sizeof(struct Source));
    double * params_hz = calloc(NP, sizeof(double));
    alloc_source(source, data->N, data->Nchannel, data->NP);
    map_array_to_params(source, params, data->T);
    map_params_to_array(source, params_hz, 1.0); 
    for(size_t i=0; i<NP; i++)
        gsl_vector_set(x,i,params_hz[i]);

    //map parameters to R
    double xmin,xmax,xn,yn, logJ = 0;
    for(size_t n=0; n<NP; n++)
    {
        xmin = gsl_matrix_get(modes[0]->minmax,n,0);
        xmax = gsl_matrix_get(modes[0]->minmax,n,1);
        xn = gsl_vector_get(x,n);
        if(xn < xmin || xn >= xmax)
        {            
            //clean up
            gsl_vector_free(x);
            free_source(source);
            return -INFINITY;
        }
        yn = logit(xn,xmin,xmax);
        gsl_vector_set(x,n,yn);
        
        //Jacobian
        logJ -= log(dsigmoid(yn, xmin, xmax));
    }
    
    //sum over modes
    double P=0.0;
    for(size_t k=0; k<NMODES; k++)
        P += modes[k]->p*multivariate_gaussian(x,modes[k]);
    
    //clean up
    gsl_vector_free(x);
    free_source(source);
    
    return log(P) + logJ;
}

double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, struct Source *source)
{
    double logP=0.0;
    double *params = source->params;
    double **uniform_prior = model->prior;
    
    //guard against nan's, but do so loudly
    if(check_range(params, uniform_prior, model->NP)) return -INFINITY;

    //update from existing runs prior
    if(flags->update)
    {
        logP = evaluate_gmm_prior(data, prior->gmm, params);
    }
    //blind search prior
    else
    {
        if(flags->volumePrior) {
            // Sky location and amplitude handled together
            logP += evaluate_volume_prior(prior, params);
            logP += evaluate_dfdtastro_prior(prior, data, model, params);
        } else {
            //sky location prior
            logP += evaluate_sky_location_prior(params, uniform_prior, model->logPriorVolume, flags->galaxyPrior, prior->skyhist, prior->dcostheta, prior->dphi, prior->nphi);

            //amplitude prior
            if(flags->snrPrior)
            {
                logP += evaluate_snr_prior(data, model, params);
            }
            else if(is_param(AMP))
            {
                if(params[AMP]<uniform_prior[AMP][0] || params[AMP]>uniform_prior[AMP][1]) return -INFINITY;
                logP -= model->logPriorVolume[AMP];
            }
        }
        
        //everything else uses simple uniform priors
        logP += evaluate_uniform_priors(params, uniform_prior, model->logPriorVolume, model->NP);
    }
    
    return logP;
}

double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume, int NP)
{
    if(check_range(params, uniform_prior, NP)) return -INFINITY;
    
    double logP = 0.0;
    //frequency bin (uniform)
    //TODO: is frequency logPriorVolume up to date?
    logP -= log(uniform_prior[F0][1]-uniform_prior[F0][0]);
    
    //cosine inclination
    logP -= logPriorVolume[COSI];
    
    //polarization
    logP -= logPriorVolume[PSI];
    
    //phase
    logP -= logPriorVolume[PHI0];
    
    //fdot (bins/Tobs)
    if(is_param(DFDT)) logP -= logPriorVolume[DFDT];
    
    //fddot
    if(is_param(D2FDT2)) logP -= logPriorVolume[D2FDT2];

    // chirp mass
    if(is_param(MC)) logP -= logPriorVolume[MC];
    
    return logP;
}

// Returns the log probability density at the location we drew
double evaluate_volume_prior(struct Prior *prior, double *params) {
    int index = sky_distance_to_sphere_index(prior, params[COSTHETA], params[PHI], params[DIST]);
    return prior->spherehist[index];
}

double evaluate_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, enum SkyPriorMode galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi)
{
    
    double logP = 0.0;
    if(galaxyFlag==galaxyPrior)
    {
        if(params[COSTHETA]<uniform_prior[COSTHETA][0] || params[COSTHETA]>uniform_prior[COSTHETA][1]) return -INFINITY;
        
        if(uniform_prior[PHI][0] > 0.0 || uniform_prior[PHI][1] < PI2)
        {
            //rejection sample on reduced prior range
            if(params[PHI]<uniform_prior[PHI][0] || params[PHI]>uniform_prior[PHI][1]) return -INFINITY;
        }
        else
        {
            //periodic boundary conditions for full range
            while(params[PHI] < 0  ) params[PHI] += PI2;
            while(params[PHI] > PI2) params[PHI] -= PI2;
        }
        //map costheta and phi to index of skyhist array
        int i = (int)floor((params[COSTHETA]-uniform_prior[COSTHETA][0])/dcostheta);
        int j = (int)floor((params[PHI]-uniform_prior[PHI][0])/dphi);
        
        int k = i*nphi + j;
        
        logP += skyhist[k];
        
        //    FILE *fptr = fopen("prior.dat","a");
        //    fprintf(fptr,"%i %i %i %g\n",i,j,k,prior->skyhist[k]);
        //    fclose(fptr);
    }
    else if(galaxyFlag==uniformPrior)
    {
        //colatitude (reflective)
        if(params[COSTHETA]<uniform_prior[COSTHETA][0] || params[COSTHETA]>uniform_prior[COSTHETA][1]) return -INFINITY;
        else logP -= logPriorVolume[COSTHETA];
        
        //longitude (periodic)
        if(uniform_prior[PHI][0] > 0.0 || uniform_prior[PHI][1] < PI2)
        {
            //rejection sample on reduced prior range
            if(params[PHI]<uniform_prior[PHI][0] || params[PHI]>uniform_prior[PHI][1]) return -INFINITY;
        }
        else
        {
            if(params[PHI]<uniform_prior[PHI][0] || params[PHI]>=uniform_prior[PHI][1])
            {
                params[PHI] = atan2(sin(params[PHI]),cos(params[PHI]));
                if(params[PHI] < 0.0) params[PHI] += PI2;
            }
        }
        logP -= logPriorVolume[PHI];
    }
    return logP;
}



double evaluate_snr_prior(struct Data *data, struct Model *model, double *params)
{
    // We should only ever call this when amplitude is a parameter in the search space
    assert(is_param(AMP));

    //check that frequency is in range
    int n = (int)floor(params[F0] - model->prior[F0][0]);
    if(n<0 || n>=data->N) return -INFINITY;
    
    //calculate noise model estimate
    double sf = data->sine_f_on_fstar;
    double sn = model->noise[0]->SnA[n]*model->noise[0]->etaA;
    
    //extra factors from TDI convention used for fractional-frequency data
    if(strcmp("frequency",data->format) == 0 || strcmp("sangria",data->format) == 0)
        sf *= asin(data->sine_f_on_fstar);

    double snr = analytic_snr(params[AMP],sn,sf,data->sqT);
    
    return log(snr_prior(snr));
}



double evaluate_calibration_prior(struct Data *data, struct Model *model)
{
    
    double dA,dphi;
    double logP = 0.0;
    
    //apply calibration error to full signal model
    //loop over time segments
    for(int m=0; m<model->NT; m++)
    {
        switch(data->Nchannel)
        {
            case 1:
                
                //amplitude
                dA   = model->calibration[m]->dampX;
                logP += log(gsl_ran_gaussian_pdf(dA,CAL_SIGMA_AMP));
                
                //phase
                dphi = model->calibration[m]->dphiX;
                logP += log(gsl_ran_gaussian_pdf(dphi,CAL_SIGMA_PHASE));
                
                break;
            case 2:
                
                //amplitude
                dA   = model->calibration[m]->dampA;
                logP += log(gsl_ran_gaussian_pdf(dA,CAL_SIGMA_AMP));
                
                //phase
                dphi = model->calibration[m]->dphiA;
                logP += log(gsl_ran_gaussian_pdf(dphi,CAL_SIGMA_PHASE));
                
                //amplitude
                dA   = model->calibration[m]->dampE;
                logP += log(gsl_ran_gaussian_pdf(dA,CAL_SIGMA_AMP));
                
                //phase
                dphi = model->calibration[m]->dphiE;
                logP += log(gsl_ran_gaussian_pdf(dphi,CAL_SIGMA_PHASE));
                
            default:
                break;
        }//end switch
    }//end loop over segments
    
    return logP;
}



