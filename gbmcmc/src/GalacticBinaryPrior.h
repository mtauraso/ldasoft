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
#include <math.h>

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
    //xcxc todo:
    // Prior ought contain more (not sure how much more) flag config information because it is passed about 
    // and packed/unpacked when creating the prior proposal.
    // It ought be that utility functions recieving a prior should not need to know --galaxy-flag or --volume-prior
    // from some other source. Information about prior configuration in general ought walk with the prior
    // Immediately this means --galaxy-prior/--volume-prior/neither state needs to be here as part of --volume-prior
    // implementation

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
 \brief Enumeration used for packing/unpacking of the prior into the proposal struct

 See setup_prior_proposal and unpack_prior_proposal for pack/unpack logic.
 This encodes the settings of flags->galaxyPrior and flags->volumePrior.
*/
enum SkyPriorMode { uniformPrior = 0, galaxyPrior = 1,  volumePrior = 2};

/**
 \brief Generate a sky prior mode enum by looking at runtime Flags
*/
enum SkyPriorMode get_sky_prior_mode(struct Flags *flags);

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
double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, struct Source *source);


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
double evaluate_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, enum SkyPriorMode galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi);

double evaluate_volume_prior(struct Prior *prior, struct Source *source);

/**
 \brief Computes uniform prior for parameters \f$\vec x\f$

 \f$ p(\vec x) = \prod_i \frac{1}{\Delta x_i} \f$
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p({\vec x})\f$
 */
double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume, int NP);
double evaluate_gmm_prior(struct Data *data, struct GMM *gmm, double *params);


/*
 These are here because they are called in inner loops, and across source files
 They ought be inlined by the compiler. In C this requires 'static inline' storage
 specifiers to ensure correct behavior between GNU C and ISO C, and that the full
 implementation be present in every translation unit (source file) the compiler processes
*/
/**
 \brief Rotate earth-origin galactic XYZ coordinates to ecliptic XYZ coordinates
*/
static inline void rotate_galtoeclip(double *xg, double *xe)
{
    xe[0] = -0.05487556043*xg[0] + 0.4941094278*xg[1] - 0.8676661492*xg[2];

    xe[1] = -0.99382137890*xg[0] - 0.1109907351*xg[1] - 0.00035159077*xg[2];

    xe[2] = -0.09647662818*xg[0] + 0.8622858751*xg[1] + 0.4971471918*xg[2];
}

/**
 \brief Rotate ecliptic XYZ coordinates to earth-origin galactic XYZ coordinates
*/
static inline void rotate_ecliptogal(double *xg, double *xe)
{
    xe[0] = -0.05487556043*xg[0] - 0.99382137890*xg[1] - 0.09647662818*xg[2];

    xe[1] =  0.4941094278*xg[0]  - 0.1109907351*xg[1]  + 0.8622858751*xg[2];

    xe[2] = -0.8676661492*xg[0]  - 0.00035159077*xg[1] + 0.4971471918*xg[2];
}

/**
 \brief Convert galactocentric XYZ to a sky position and distance.
*/
static inline void galactocentric_to_sky_distance(/*in*/ double x[3], /*out*/ double *phi, /*out*/ double *theta, /*out*/ double *r_ec)
{
    double xe[3], xg[3];

    xg[0] = x[0] - GALAXY_RGC;   // solar barycenter is offset from galactic center along x-axis (by convention)
    xg[1] = x[1];
    xg[2] = x[2];

    /* Rotate from galactic to ecliptic */
    rotate_galtoeclip(xg, xe);

    *r_ec = sqrt(xe[0]*xe[0]+xe[1]*xe[1]+xe[2]*xe[2]);

    *theta = M_PI/2.0-acos(xe[2]/(*r_ec));

    *phi = atan2(xe[1],xe[0]);

    if(*phi<0.0) *phi += 2.0*M_PI;
}

/**
 \brief Convert sky position and distance to a galactocentric XYZ
*/
static inline void sky_distance_to_galactocentric(/*out*/ double x[3], /*in*/ double phi, /*in*/ double theta, /*in*/ double r_ec)
{
    double xe[3], xg[3];

    // Shift theta to be azimuthal angle to z axis rather than an elevation from xy plane
    theta = M_PI/2.0 - theta;

    // Convert to ecliptic xyz
    xe[0] = r_ec * cos(phi) * sin(theta);
    xe[1] = r_ec * sin(phi) * sin(theta);
    xe[2] = r_ec * cos(theta);

    // Rotate to galactic
    rotate_ecliptogal(xe, xg);

    // Solar barycenter is offset from galactic center along x axis
    x[0] = xg[0] + GALAXY_RGC;
    x[1] = xg[1];
    x[2] = xg[2];
}

/**
\brief Uniform draw of sample location within the galactic bounding box
 */
static inline void _generate_uniform_galaxy_sample(/*out*/ double * proposed_sample, gsl_rng *r)
{
    proposed_sample[0] = (GALAXY_BB_X*0.5)*(-1.0+2.0*gsl_rng_uniform(r));
    proposed_sample[1] = (GALAXY_BB_Y*0.5)*(-1.0+2.0*gsl_rng_uniform(r));
    proposed_sample[2] = (GALAXY_BB_Z*0.5)*(-1.0+2.0*gsl_rng_uniform(r));
}

/**
\brief Uniform draw of sample location within the galactic bounding box, and return the log of the volume drawn from.
 */
static inline double generate_uniform_galaxy_sample(/*out*/ double * proposed_sample, gsl_rng *r)
{
    _generate_uniform_galaxy_sample(proposed_sample, r);
    return  log(GALAXY_BB_X) + log(GALAXY_BB_Y) + log(GALAXY_BB_Z);
}

#endif /* GalacticBinaryPrior_h */
