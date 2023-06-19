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
 @file GalacticBinaryPrior.h
 \brief Functions supporting prior distributions.
 
 Includes functions for createing and evaluating model priors.
 */


#ifndef GalacticBinaryPrior_h
#define GalacticBinaryPrior_h

///@name Galaxy Prior
///@{
#define GALAXY_RGC 7.2 //!<distance from solar BC to GC (kpc)
#define GALAXY_A  0.25 //!<bulge fraction
#define GALAXY_Rb 0.8  //!< bulge radius (kpc)
#define GALAXY_Rd 2.5  //!< disk radius (kpc)
#define GALAXY_Zd 0.4  //!< disk height (kpc)

#define GALAXY_BB_X (40.0*GALAXY_Rd) //!< Dimension of galactic bounding box in galactocentric x (kpc) used in generating sky location and 3d galaxy prior.
#define GALAXY_BB_Y (40.0*GALAXY_Rd) //!< Dimension of galactic bounding box in galactocentric y (kpc) used in generating sky location and 3d galaxy prior.
#define GALAXY_BB_Z (80.0*GALAXY_Zd) //!< Dimension of galactic bounding box in galactocentric z (kpc) used in generating sky location and 3d galaxy prior.
///@}

///@name Calibration prior
///@{
#define CAL_SIGMA_PHASE 0.35 //!< 1-\f$\sigma\f$ phase error (rad) \f${\sim}30^\circ\f$
#define CAL_SIGMA_AMP 0.20 //!< 1-\f$\sigma\f$ fractional amplitude error
///@}

/*!
 \brief Prototype structure for prior distributions.
 
 Generic data structure for holding all information needed by prior distributions.
 Structure contains parameters for different supported priors and various book-keeping scalars, vectors, and matrices to
 hold needed metadata.

*/
struct Prior
{    
    ///@name Uniform prior
    ///@{
    double **prior; //!<upper and lower bounds for uniform priors \f$ [\theta_{\rm min},\theta_{\rm max}]\f$
    double logPriorVolume; //!<prior volume \f$ -\sum \log(\theta_{\rm max}-\theta_{\rm min})\f$
    ///@}

    ///@name Sky Location prior
    ///@{
    double *skyhist; //!<2D histogram of prior density on sky
    double dcostheta; //!<size of `skyhist` bins in \f$\cos\theta\f$ direction
    double dphi; //!<size of `skyhist` bins in \f$\phi\f$ direction
    double skymaxp; //!<max prior density of `skyhist`
    int ncostheta; //!<number of `skyhist` bins in \f$\cos\theta\f$ direction
    int nphi; //!<number of `skyhist` bins in \f$\phi\f$ direction
    ///@}

    ///@name 3D galaxy prior
    ///@{
    double * volhist; //!< 3D histogram of prior density in galactocentric coordinates (x,y,z)
    double volmaxp; //!< Max prior density of volhist
    int nx; //!< number of bins in x
    int ny; //!< number of bins in y
    int nz; //!< number of bins in z
    double dx; //!< size of a bin in kpc x direction
    double dy; //!< size of a bin in kpc y direction
    double dz; //!< size of a bin in kpc z direction
    ///@}
    
    ///@name workspace
    ///@{
    double *vector;  //!<utility 1D array for prior metadata
    double **matrix; //!<utility 2D array for prior metadata
    double ***tensor;//!<utility 3D array for prior metadata
    ///@}

    /// Gaussian Mixture Model prior
    struct GMM *gmm;
};
/**
 \brief Checks that parameters \f$\vec x\f$ are within prior volume \f$V\f$.
 
 @returns 0 if \f$\vec x \in V\f$
 @returns 1 else
 */
int check_range(double *params, double **uniform_prior, int NP);

/**
\brief Set up sky location prior assuming galactic distribution
 
 Uses axially symmetric disk and bulge model of galaxy,
 randomly distributes points within the assumed distribution,
 and bins the points based on their sky location from Earth
 to use as a prior.
 */
void set_galaxy_prior(struct Flags *flags, struct Prior *prior);

/**
\brief Sets Gaussian Mixture Model prior for source model
 */
void set_gmm_prior(struct Flags *flags, struct Data *data, struct Prior *prior);

/**
 \brief Sets Uniform prior for source model
 */
void set_uniform_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose);

/**
 \brief Allocates a prior structure that can be used with set_ functions
*/
struct Prior * alloc_prior();

/**
 \brief Frees any previously allocated prior structure
*/
void free_prior(struct Prior *prior);

/**
 \brief Computes joint prior for input parameters `params`
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p(\vec x)\f$
 */
double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params);


/**
 \brief Computes prior amplitude parameter `params`
 
 Uses analytic approximation to signal to noise ratio
 and computes prior via snr_prior().
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p({\rm SNR})\f$
 */
double evaluate_snr_prior(struct Data *data, struct Model *model, double *params);

/**
 \brief Computes prior for sky location parameters \f$\vec\Omega\f$
 
 Depending on Flags::galaxyFlag, evaluates either uniform or galaxy prior
 for sky location parameters
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p({\vec\Omega})\f$
 */
double evalaute_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, int galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi);

/**
 \brief Computes uniform prior for parameters \f$\vec x\f$

 \f$ p(\vec x) = \prod_i \frac{1}{\Delta x_i} \f$
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p({\vec x})\f$
 */
double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume, int NP);
double evaluate_gmm_prior(struct Data *data, struct GMM *gmm, double *params);



#endif /* GalacticBinaryPrior_h */
