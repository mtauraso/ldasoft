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

/**
 @file GalacticBinaryModel.h
 \brief Functions for defining, manipulating, and evaluating Galactic Binary model
 */


#ifndef GalacticBinaryModel_h
#define GalacticBinaryModel_h

#include <stdio.h>

/**
\brief Create galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and creates meta template of all sources in the model.
 */
void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
\brief F-statistic maximization of galactic binary parameters

 Wrapper of GalacticBinaryFstatistic.c functions for maximizing waveform
 over \f$(\mathcal{A},\cos\iota,\psi,\varphi_0)\f$ by filtering on
 original data.
 
 \todo test filtering on residuals
 */
void maximize_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
 \brief Create LISA instrument noise model
 
 Computes \f$S_n(f)\f$ from Model::noise.
 */
void generate_noise_model(struct Data *data, struct Model *model);

/**
 \brief Create LISA instrument calibration model
 
 Computes phase and amplitude corrections from Model::calibration parameters.
 */
void generate_calibration_model(struct Data *data, struct Model *model);

/**
\brief Apply amplitude and phase corrections.

Computes new LISA instrument response Model::tdi after applying
 amplitude and phase corrections from calibration parameters.
 */
void apply_calibration_model(struct Data *data, struct Model *model);

/**
 \brief Compute argument of Gaussian likelihood
 
 Computes residual of data and meta-template from Model and noise weighted inner product.
 @return \f$ -\frac{1}{2}(d-h|d-h) \f$
 */
double gaussian_log_likelihood(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for constant noise level
 
 For noise models that are (approximately) a constant over the band
 of \f$N\f$ Fourier bins, parameterized by a multiplyer \f$\eta\f$
 
 @return \f$ \sum_{\rm TDI} N\log\eta_{\rm TDI} \f$
 */
double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for arbitrary noise level
 
 For noise models that are free to vary over the band
 of \f$N\f$ Fourier bins
 
 @return \f$ \sum_{\rm TDI} \sum_f \log S_{n,{\rm TDI}}(f) \f$
 */

double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model);

/**
 \brief Check for increase in maximum log likelihood
 */
int update_max_log_likelihood(struct Model **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Converts physical UCB parameters to array expected by GalacticBinaryWaveform.c
 */
static inline void map_params_to_array(struct Source *source, double *params, double T);

/**
 \brief Converts array expected by GalacticBinaryWaveform.c to
 physical UCB parameters
 */
static inline void map_array_to_params(struct Source *source, double *params, double T);

/**
 \brief Loads spacecraft orbits, either analytically or from tabulated ephemerides.
 */
void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Allocates and initializes Chain structure and prepares output files.
 */
void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, const char *mode);

/** @name Allocate memory for structures */
///@{
void alloc_data(struct Data *data, struct Flags *flags);
void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel, int NP, int NT);
void alloc_noise(struct Noise *noise, int NFFT);
void realloc_noise(struct Noise *noise, int NFFT);
void alloc_source(struct Source *source, int NFFT, int Nchannel, int NP);
void alloc_calibration(struct Calibration *calibration);
///@}

/**
 \brief Shallow copy of Data structure
 */
void copy_data(struct Data *origin, struct Data *copy);

/**
 \brief Very Shallow copy of Source. 
        Only touches params array, NP, and named parameter members which can be calculated from params.
 */
void copy_source_params_only(struct Source *origin, struct Source *copy, double T);

/** @name Deep copy structure contents */
///@{
void copy_source(struct Source *origin, struct Source *copy);
void copy_model(struct Model *origin, struct Model *copy);
void copy_noise(struct Noise *origin, struct Noise *copy);
void copy_model_lite(struct Model *origin, struct Model *copy);
void copy_calibration(struct Calibration *origin, struct Calibration *copy);
///@}

/** @name Free memory for structures */
///@{
void free_noise(struct Noise *noise);
void free_model(struct Model *model);
void free_source(struct Source *source);
void free_chain(struct Chain *chain, struct Flags *flags);
void free_calibration(struct Calibration *calibration);
///@}

/**
 \brief Deep comparison of Model contents
 
 @return 0 if models are identical
 @return 1 if models are different
 */
int compare_model(struct Model *a, struct Model *b);

// Definitions in the header file because these are used all over the place for
// parameter-layout agnostic unit changes in frequency space, or simply picking out
// mandatory parameters from params array. 
//
// This is done to allow the compiler to optimize out unnecessary math on variables
// that will go unused in the scope of the caller.

#include <assert.h>  // for assert

#include "GalacticBinary.h" // for struct Source
#include "GalacticBinaryWaveform.h" // for galactic_binary_Amp and galactic_binary_fdot

static inline void map_array_to_params(struct Source *source, double *params, double T)
{
    source->f0       = params[F0]/T;
    source->costheta = params[COSTHETA];
    source->phi      = params[PHI];

    // All legacy params are needed for file I/O because sources and catalog files must 
    // reflect the parameters of the gravitational wave, regardless of what 
    // parameterization is in use in memory.
    //
    // The legacy params are F0, COSTHETA, PHI, AMP, COSI, PSI, PHI0, and DFDT
    // Therefore we must always compute them. This assert is paranoia to that effect
    // which matches if statements below.
    //
    // We do these as asserts to prevent unnecessary branches in release code, because
    // this function is called a lot.
    //
    // Assure we calculate amplitude
    assert(is_param(AMP) || (is_param(DIST) && is_param(MC)));
    // Assure we calculate dfdt
    assert(is_param(DFDT) || (is_param(DFDTASTRO) && is_param(MC)));

    if(is_param(AMP))
        source->amp      = exp(params[AMP]);
        // Note that when in amplitude only mode there is no
        // unambiguous way to calculate Mc and D, without assuming a
        // detached binary.
    if(is_param(DIST) && is_param(MC)) {
        source->amp = galactic_binary_Amp(params[MC], source->f0, params[DIST]);
        source->Mc  = params[MC];
        source->D   = params[DIST];
    }

    source->cosi     = params[COSI];
    source->psi      = params[PSI];
    source->phi0     = params[PHI0];

    source->dfdt   = 0.0;  // This line is slightly paranoid given the assert about DFDT above
    source->d2fdt2 = 0.0;
    if(is_param(DFDT))
        source->dfdt   = params[DFDT]/(T*T);
    if(is_param(DFDTASTRO) && is_param(MC)) { 
        source->dfdtastro = params[DFDTASTRO]/(T*T);
        source->dfdt   =  source->dfdtastro + galactic_binary_fdot(params[MC], source->f0);
    }
    if(is_param(D2FDT2))
        source->d2fdt2 = params[D2FDT2]/(T*T*T);
}

static inline void map_params_to_array(struct Source *source, double *params, double T)
{
    params[F0] = source->f0*T;
    params[COSTHETA] = source->costheta;
    params[PHI] = source->phi;
    if(is_param(AMP))
        params[AMP] = log(source->amp);
    params[COSI] = source->cosi;
    params[PSI] = source->psi;
    params[PHI0] = source->phi0;
    if(is_param(DFDT))
        params[DFDT] = source->dfdt*T*T;
    if(is_param(DIST))
        params[DIST] = source->D;
    if(is_param(MC))
        params[MC] = source->Mc;
    if(is_param(DFDTASTRO) && is_param(MC))
        params[DFDTASTRO] = (source->dfdt - galactic_binary_fdot(source->Mc, source->f0))*T*T;
    if(is_param(D2FDT2))
        params[D2FDT2] = source->d2fdt2*T*T*T;
}

#endif /* GalacticBinaryModel_h */
