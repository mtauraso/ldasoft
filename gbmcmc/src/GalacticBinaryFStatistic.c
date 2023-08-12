/*
 *  Copyright (C) 2019 Travis Robson, Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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


//
//  GalacticBinaryFStatistic.c
//  
//
//  Created on 7/21/17 by
//    Robson, Travis (Montana State Univ.)
//    Cornish, Neil (Montana State Univ.)
//    Littenberg, Tyson B. (MSFC-ZP12)
//
//
#include <math.h>

#include <LISA.h>

#include "GalacticBinary.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryFStatistic.h"


void initialize_XLS(long M, double *XLS, double *AA, double *EE)
{
    for (int i=0; i<2*M; i++) // initialize arrays to zeros
    {
        XLS[i] = 0.0;
        AA[i]  = 0.0;
        EE[i]  = 0.0;
    }
}

void init_A_filters(struct Orbit *orbit, struct Data *data, struct Filter *F_filter)
{
    long M = F_filter->M_filter;
    
    F_filter->A1_fX = calloc((2*M),sizeof(double));
    F_filter->A2_fX = calloc((2*M),sizeof(double));
    F_filter->A3_fX = calloc((2*M),sizeof(double));
    F_filter->A4_fX = calloc((2*M),sizeof(double));
    
    F_filter->A1_fA = calloc((2*M),sizeof(double));
    F_filter->A2_fA = calloc((2*M),sizeof(double));
    F_filter->A3_fA = calloc((2*M),sizeof(double));
    F_filter->A4_fA = calloc((2*M),sizeof(double));
    
    F_filter->A1_fE = calloc((2*M),sizeof(double));
    F_filter->A2_fE = calloc((2*M),sizeof(double));
    F_filter->A3_fE = calloc((2*M),sizeof(double));
    F_filter->A4_fE = calloc((2*M),sizeof(double));
    
    initialize_XLS(M, F_filter->A1_fX, F_filter->A1_fA, F_filter->A1_fE);
    initialize_XLS(M, F_filter->A2_fX, F_filter->A2_fA, F_filter->A2_fE);
    initialize_XLS(M, F_filter->A3_fX, F_filter->A3_fA, F_filter->A3_fE);
    initialize_XLS(M, F_filter->A4_fX, F_filter->A4_fA, F_filter->A4_fE);
    
    //   get_filters(orbit, data, 1, F_filter);
    //   get_filters(orbit, data, 2, F_filter);
    //   get_filters(orbit, data, 2, F_filter);
    //   get_filters(orbit, data, 2, F_filter);
    
    // Make use of a phase shift to quickly generate other filters
    get_filters(orbit, data, 3, F_filter);    
    // copy  F_filter->A3_fX into F_filter->A1_fX
    for (int i=0; i<M; i++)
    {
        F_filter->A1_fX[2*i+1]   = -F_filter->A3_fX[2*i];
        F_filter->A1_fX[2*i] =  F_filter->A3_fX[2*i+1];
        
        F_filter->A1_fA[2*i+1]   = -F_filter->A3_fA[2*i];
        F_filter->A1_fA[2*i] =  F_filter->A3_fA[2*i+1];
        
        F_filter->A1_fE[2*i+1]   = -F_filter->A3_fE[2*i];
        F_filter->A1_fE[2*i] =  F_filter->A3_fE[2*i+1];
    }
    
    get_filters(orbit, data, 4, F_filter);
    for (int i=0; i<M; i++)
    {
        F_filter->A2_fX[2*i+1]    = -F_filter->A4_fX[2*i];
        F_filter->A2_fX[2*i] =  F_filter->A4_fX[2*i+1];
        
        F_filter->A2_fA[2*i+1]    = -F_filter->A4_fA[2*i];
        F_filter->A2_fA[2*i] =  F_filter->A4_fA[2*i+1];
        
        F_filter->A2_fE[2*i+1]    = -F_filter->A4_fE[2*i];
        F_filter->A2_fE[2*i] =  F_filter->A4_fE[2*i+1];
    }
    
    //////////////////////
}

void init_M_matrix(struct Filter *F_filter, struct Data *data)
{
    long i,j;
    
    F_filter->M_inv_X  = malloc(4*sizeof(double *));//dmatrix(0,3,0,3);
    F_filter->M_inv_AE = malloc(4*sizeof(double *));//dmatrix(0,3,0,3);
    for(i=0; i<4; i++)
    {
        F_filter->M_inv_X[i]  = calloc(4,sizeof(double));//dmatrix(0,3,0,3);
        F_filter->M_inv_AE[i] = calloc(4,sizeof(double));//dmatrix(0,3,0,3);
        
        for (j=0;j<4;j++)
        {
            F_filter->M_inv_X[i][j]  = 0.0;
            F_filter->M_inv_AE[i][j] = 0.0;
        }
    }
    
    get_M(F_filter, F_filter->M_inv_X, F_filter->M_inv_AE, data);
    
}

void free_Filter(struct Filter *F_filter)
{
    free(F_filter->A1_fX);
    free(F_filter->A2_fX);
    free(F_filter->A3_fX);
    free(F_filter->A4_fX);
    
    free(F_filter->A1_fA);
    free(F_filter->A2_fA);
    free(F_filter->A3_fA);
    free(F_filter->A4_fA);
    
    free(F_filter->A1_fE);
    free(F_filter->A2_fE);
    free(F_filter->A3_fE);
    free(F_filter->A4_fE);
    
    for(int i=0; i<4; i++)
    {
        free(F_filter->M_inv_X[i]);
        free(F_filter->M_inv_AE[i]);
    }
    free(F_filter->M_inv_X);
    free(F_filter->M_inv_AE);
    free(F_filter);
}

void get_filters(struct Orbit *orbit, struct Data *data, int filter_id, struct Filter *F_filter)
{
    struct Source source;
    long M_filter;
    double *X, *A, *E;
    
    M_filter = F_filter->M_filter;
    
    source.f0 = F_filter->f0;				 // f0 
    source.costheta = cos(F_filter->theta); // theta
    source.phi = F_filter->phi;    	     // phi
    
    source.dfdt = F_filter->fdot;				 // fdot
    source.d2fdt2 = F_filter->fddot;			 // fddot
    
    // Now need to calculate the Filters A^{i} (Cornish & Crowder '05)
    // A^{1} -> iota = pi/2, psi = 0,    A = 2, phase = 0
    // A^{2} -> iota = pi/2, psi = pi/4, A = 2, phase = pi
    // A^{3} -> iota = pi/2, psi = 0,    A = 2, phase = 3*pi/2
    // A^{4} -> iota = pi/2, psi = pi/4, A = 2, phase = pi/2

    // Same values across all filters
    source.amp = 2.0;
    source.cosi = cos(PIon2);
    
    // determine which parameters to pass into galactic_binary
    if (filter_id == 1)
    {    
        source.psi = 0.0;
        source.phi0 = 0.0;
        X = F_filter->A1_fX;
        A = F_filter->A1_fA;
        E = F_filter->A1_fE;
    } else if (filter_id == 2){
        source.psi = PIon4;
        source.phi0 = M_PI;
        X = F_filter->A2_fX;
        A = F_filter->A2_fA;
        E = F_filter->A2_fE;
    } else if (filter_id == 3){
        source.psi = 0.0;
        source.phi0 = 3.0*PIon2;
        X = F_filter->A3_fX;
        A = F_filter->A3_fA;
        E = F_filter->A3_fE;
    } else {
        source.psi = PIon4;
        source.phi0 = PIon2;
        X = F_filter->A4_fX;
        A = F_filter->A4_fA;
        E = F_filter->A4_fE;
    }
    galactic_binary_from_source(orbit, data->format, data->T, data->t0[0], &source, data->NP, X, A, E, M_filter, 2);
}

void get_N(struct Data *data, struct Filter *F_filter)
{
    long     i,k;
    long     M_filter;
    long q;
    double *XfLS, *AALS, *EELS;
    
    /////////
    //
    // Calculate N^{i} (Cornish & Crowder '05)
    // N^{i} = (s|A^{i})
    //
    /////////
    XfLS = data->tdi[FIXME]->X;
    AALS = data->tdi[FIXME]->A;
    EELS = data->tdi[FIXME]->E;
    
    q  = F_filter->q - data->qmin;
    
    M_filter = F_filter->M_filter;
    
    F_filter->N1_X = 0.; F_filter->N1_AE = 0.;
    F_filter->N2_X = 0.; F_filter->N2_AE = 0.;
    F_filter->N3_X = 0.; F_filter->N3_AE = 0.;
    F_filter->N4_X = 0.; F_filter->N4_AE = 0.;
    
    for (i=0; i<M_filter; i++)
    {
        k = (q + i - M_filter/2);
        
        if(k>0 && k<data->N)
        {
            
            F_filter->N1_X  += (XfLS[2*k]*F_filter->A1_fX[2*i] + XfLS[2*k+1]*F_filter->A1_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            F_filter->N2_X  += (XfLS[2*k]*F_filter->A2_fX[2*i] + XfLS[2*k+1]*F_filter->A2_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            F_filter->N3_X  += (XfLS[2*k]*F_filter->A3_fX[2*i] + XfLS[2*k+1]*F_filter->A3_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            F_filter->N4_X  += (XfLS[2*k]*F_filter->A4_fX[2*i] + XfLS[2*k+1]*F_filter->A4_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            
            F_filter->N1_AE += (AALS[2*k]*F_filter->A1_fA[2*i] + AALS[2*k+1]*F_filter->A1_fA[2*i+1]
                                +EELS[2*k]*F_filter->A1_fE[2*i] + EELS[2*k+1]*F_filter->A1_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            F_filter->N2_AE += (AALS[2*k]*F_filter->A2_fA[2*i] + AALS[2*k+1]*F_filter->A2_fA[2*i+1]
                                +EELS[2*k]*F_filter->A2_fE[2*i] + EELS[2*k+1]*F_filter->A2_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            F_filter->N3_AE += (AALS[2*k]*F_filter->A3_fA[2*i] + AALS[2*k+1]*F_filter->A3_fA[2*i+1]
                                +EELS[2*k]*F_filter->A3_fE[2*i] + EELS[2*k+1]*F_filter->A3_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            F_filter->N4_AE += (AALS[2*k]*F_filter->A4_fA[2*i] + AALS[2*k+1]*F_filter->A4_fA[2*i+1]
                                +EELS[2*k]*F_filter->A4_fE[2*i] + EELS[2*k+1]*F_filter->A4_fE[2*i+1])/data->noise[FIXME]->SnA[k];
        }
    }
    
    F_filter->N1_X  *= 4.0;
    F_filter->N2_X  *= 4.0;
    F_filter->N3_X  *= 4.0;
    F_filter->N4_X  *= 4.0;
    
    F_filter->N1_AE *= 4.0;
    F_filter->N2_AE *= 4.0;
    F_filter->N3_AE *= 4.0;
    F_filter->N4_AE *= 4.0;
    
}

void get_M(struct Filter *F_filter, double **M_inv_X, double **M_inv_AE, struct Data *data)
{
    int i,j,k;
    long M;
    long q, M_filter;
    
    M = F_filter->M_filter;
    q = F_filter->q-data->qmin;
    M_filter = F_filter->M_filter;
    
    for (i=0; i<M; i++)
    {
        k = (q + i - M_filter/2);
        
        if(k>0 && k<data->N)
        {
            M_inv_X[0][0] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A1_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A1_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            M_inv_X[0][1] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A2_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A2_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            M_inv_X[0][2] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A3_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A3_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            M_inv_X[0][3] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            
            M_inv_X[1][1] += 4.0*(F_filter->A2_fX[2*i]  *F_filter->A2_fX[2*i]
                                  +F_filter->A2_fX[2*i+1]*F_filter->A2_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            M_inv_X[1][2] += 4.0*(F_filter->A2_fX[2*i]  *F_filter->A3_fX[2*i]
                                  +F_filter->A2_fX[2*i+1]*F_filter->A3_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            M_inv_X[1][3] += 4.0*(F_filter->A2_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A2_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            
            M_inv_X[2][2] += 4.0*(F_filter->A3_fX[2*i]  *F_filter->A3_fX[2*i]
                                  +F_filter->A3_fX[2*i+1]*F_filter->A3_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            M_inv_X[2][3] += 4.0*(F_filter->A3_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A3_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            
            M_inv_X[3][3] += 4.0*(F_filter->A4_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A4_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise[FIXME]->SnX[k];
            
            
            
            
            M_inv_AE[0][0] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A1_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A1_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A1_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A1_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            M_inv_AE[0][1] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A2_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A2_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A2_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A2_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            M_inv_AE[0][2] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A3_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A3_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A3_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A3_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            M_inv_AE[0][3] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            
            
            M_inv_AE[1][1] += 4.0*(F_filter->A2_fA[2*i]  *F_filter->A2_fA[2*i]
                                   +F_filter->A2_fA[2*i+1]*F_filter->A2_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A2_fE[2*i]  *F_filter->A2_fE[2*i]
                  +F_filter->A2_fE[2*i+1]*F_filter->A2_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            M_inv_AE[1][2] += 4.0*(F_filter->A2_fA[2*i]  *F_filter->A3_fA[2*i]
                                   +F_filter->A2_fA[2*i+1]*F_filter->A3_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A2_fE[2*i]  *F_filter->A3_fE[2*i]
                  +F_filter->A2_fE[2*i+1]*F_filter->A3_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            M_inv_AE[1][3] += 4.0*(F_filter->A2_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A2_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A2_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A2_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            
            
            M_inv_AE[2][2] += 4.0*(F_filter->A3_fA[2*i]  *F_filter->A3_fA[2*i]
                                   +F_filter->A3_fA[2*i+1]*F_filter->A3_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A3_fE[2*i]  *F_filter->A3_fE[2*i]
                  +F_filter->A3_fE[2*i+1]*F_filter->A3_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            M_inv_AE[2][3] += 4.0*(F_filter->A3_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A3_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A3_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A3_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise[FIXME]->SnA[k];
            
            
            M_inv_AE[3][3] += 4.0*(F_filter->A4_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A4_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise[FIXME]->SnA[k]
            +4.0*(F_filter->A4_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A4_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise[FIXME]->SnA[k];
        }
    }
    
    // make use of symmetry to populate the other pieces
    for (j=0; j<4; j++)
    {
        for (k=j+1; k<4; k++)
        {
            M_inv_X[k][j]  = M_inv_X[j][k];
            M_inv_AE[k][j] = M_inv_AE[j][k];
        }
    }
    
    invert_matrix(M_inv_X,4);
    invert_matrix(M_inv_AE,4);
    
}


int sgn(double v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

void get_F_params(struct Filter *F_filter)
{
    
    double Aplus_X,Across_X;
    double Aplus_AE,Across_AE;
    
    double psi_X_Fstat,  c,   A_X_Fstat,  iota_X_Fstat,  phase_X_Fstat;
    double psi_AE_Fstat, cAE, A_AE_Fstat, iota_AE_Fstat, phase_AE_Fstat;
    
    // TODO: redundant terms here, save time break it up
    double a1p4X = F_filter->a1_X + F_filter->a4_X;
    double a2m3X = F_filter->a2_X - F_filter->a3_X;
    double a1m4X = F_filter->a1_X - F_filter->a4_X;
    double a2p3X = F_filter->a2_X + F_filter->a3_X;
    
    Aplus_X  = sqrt( a1p4X*a1p4X + a2m3X*a2m3X  ) + sqrt( a1m4X*a1m4X + a2p3X*a2p3X );
    Across_X = sqrt( a1p4X*a1p4X + a2m3X*a2m3X  ) - sqrt( a1m4X*a1m4X + a2p3X*a2p3X );
    
    psi_X_Fstat   = 0.5*atan2(Aplus_X*F_filter->a4_X - Across_X*F_filter->a1_X,
                              -(Across_X*F_filter->a2_X + Aplus_X*F_filter->a3_X) );
    
    c             = (double)sgn(sin(2.0*psi_X_Fstat));
    A_X_Fstat     = 0.5*(   Aplus_X + sqrt(Aplus_X*Aplus_X - Across_X*Across_X)  );
    iota_X_Fstat  = acos(  -Across_X/(Aplus_X + sqrt(Aplus_X*Aplus_X - Across_X*Across_X))  );
    phase_X_Fstat = atan2(  c*(Aplus_X*F_filter->a4_X - Across_X*F_filter->a1_X),
                          -c*(Across_X*F_filter->a3_X + Aplus_X*F_filter->a2_X)  );
    
    double a1p4AE = F_filter->a1_AE + F_filter->a4_AE;
    double a2m3AE = F_filter->a2_AE - F_filter->a3_AE;
    double a1m4AE = F_filter->a1_AE - F_filter->a4_AE;
    double a2p3AE = F_filter->a2_AE + F_filter->a3_AE;
    
    Aplus_AE  = sqrt( a1p4AE*a1p4AE + a2m3AE*a2m3AE  ) + sqrt( a1m4AE*a1m4AE + a2p3AE*a2p3AE );
    Across_AE = sqrt( a1p4AE*a1p4AE + a2m3AE*a2m3AE  ) - sqrt( a1m4AE*a1m4AE + a2p3AE*a2p3AE );
    
    
    psi_AE_Fstat   = atan2(Aplus_AE*F_filter->a4_AE - Across_AE*F_filter->a1_AE,
                               -(Across_AE*F_filter->a2_AE + Aplus_AE*F_filter->a3_AE) );
    
    psi_AE_Fstat *= 0.5;
    
    cAE            = (double)sgn(sin(2.0*psi_AE_Fstat));
    A_AE_Fstat     = 0.5*(  Aplus_AE + sqrt(Aplus_AE*Aplus_AE - Across_AE*Across_AE)  );
    iota_AE_Fstat  = acos(  -Across_AE/(Aplus_AE + sqrt(Aplus_AE*Aplus_AE - Across_AE*Across_AE))  );
    phase_AE_Fstat = atan2(  cAE*(Aplus_AE*F_filter->a4_AE - Across_AE*F_filter->a1_AE),
                           -cAE*(Across_AE*F_filter->a3_AE + Aplus_AE*F_filter->a2_AE)  );

    // 	printf("\nX F-stat ML Values                         AE F-stat ML Values\n");
    // 	printf("-------------------------------            -------------------------------\n");
    // 	printf("psi X:   	%e               psi AE:   	%e\n",  tan(psi_X_Fstat),   tan(psi_AE_Fstat));
    // 	printf("A X:     	%e               A AE:     	%e\n",  A_X_Fstat,          A_AE_Fstat);
    // 	printf("iota X:  	%e               iota AE:  	%e \n", cos(iota_X_Fstat),  cos(iota_AE_Fstat));
    // 	printf("phase X: 	%e               phase AE: 	%e\n",  tan(phase_X_Fstat), tan(phase_AE_Fstat));
    //
    // 	printf("\nDifference                                 Difference\n");
    // 	printf("-------------------------------            -------------------------------\n");
    // 	printf("psi X:   	%e               psi AE:   	%e\n", tan(psi_X_Fstat)-tan(src->params[5]),   tan(psi_AE_Fstat)-tan(src->params[5]));
    // 	printf("A X:     	%e               A AE:     	%e\n", A_X_Fstat - src->params[3],             A_AE_Fstat - src->params[3]);
    // 	printf("iota X:  	%e               iota AE:  	%e\n", cos(iota_X_Fstat)-cos(src->params[4]),  cos(iota_AE_Fstat)-cos(src->params[4]));
    // 	printf("phase X: 	%e              phase AE: 	%e\n", tan(phase_X_Fstat)-tan(src->params[6]), tan(phase_AE_Fstat)-tan(src->params[6]));
    //
    // 	printf("\n================================================================================\n\n");
    
    F_filter->psi_X_Fstat   = psi_X_Fstat;
    F_filter->A_X_Fstat     = A_X_Fstat;
    F_filter->iota_X_Fstat  = iota_X_Fstat;
    F_filter->phase_X_Fstat = phase_X_Fstat;
    
    F_filter->psi_AE_Fstat   = psi_AE_Fstat;
    F_filter->A_AE_Fstat     = A_AE_Fstat;
    F_filter->iota_AE_Fstat  = iota_AE_Fstat;
    F_filter->phase_AE_Fstat = phase_AE_Fstat;
}

void calc_Fstat_logL(struct Filter *F_filter)
{
    F_filter->Fstat_X  = 0.5*(F_filter->M_inv_X[0][0]*F_filter->N1_X*F_filter->N1_X
                              + 2*F_filter->M_inv_X[0][1]*F_filter->N1_X*F_filter->N2_X
                              + 2*F_filter->M_inv_X[0][2]*F_filter->N1_X*F_filter->N3_X
                              + 2*F_filter->M_inv_X[0][3]*F_filter->N1_X*F_filter->N4_X
                              +   F_filter->M_inv_X[1][1]*F_filter->N2_X*F_filter->N2_X
                              + 2*F_filter->M_inv_X[1][2]*F_filter->N2_X*F_filter->N3_X
                              + 2*F_filter->M_inv_X[1][3]*F_filter->N2_X*F_filter->N4_X
                              +   F_filter->M_inv_X[2][2]*F_filter->N3_X*F_filter->N3_X
                              + 2*F_filter->M_inv_X[2][3]*F_filter->N3_X*F_filter->N4_X
                              +   F_filter->M_inv_X[3][3]*F_filter->N4_X*F_filter->N4_X);
    
    F_filter->Fstat_AE = 0.5*(F_filter->M_inv_AE[0][0]*F_filter->N1_AE*F_filter->N1_AE
                              + 2*F_filter->M_inv_AE[0][1]*F_filter->N1_AE*F_filter->N2_AE
                              + 2*F_filter->M_inv_AE[0][2]*F_filter->N1_AE*F_filter->N3_AE
                              + 2*F_filter->M_inv_AE[0][3]*F_filter->N1_AE*F_filter->N4_AE
                              +   F_filter->M_inv_AE[1][1]*F_filter->N2_AE*F_filter->N2_AE
                              + 2*F_filter->M_inv_AE[1][2]*F_filter->N2_AE*F_filter->N3_AE
                              + 2*F_filter->M_inv_AE[1][3]*F_filter->N2_AE*F_filter->N4_AE
                              +   F_filter->M_inv_AE[2][2]*F_filter->N3_AE*F_filter->N3_AE
                              + 2*F_filter->M_inv_AE[2][3]*F_filter->N3_AE*F_filter->N4_AE
                              +   F_filter->M_inv_AE[3][3]*F_filter->N4_AE*F_filter->N4_AE);
}

void calc_a_i(struct Filter *F_filter)
{
    
    F_filter->a1_X = F_filter->M_inv_X[0][0]*F_filter->N1_X + F_filter->M_inv_X[0][1]*F_filter->N2_X
    +F_filter->M_inv_X[0][2]*F_filter->N3_X + F_filter->M_inv_X[0][3]*F_filter->N4_X;
    F_filter->a2_X = F_filter->M_inv_X[1][0]*F_filter->N1_X + F_filter->M_inv_X[1][1]*F_filter->N2_X
    +F_filter->M_inv_X[1][2]*F_filter->N3_X + F_filter->M_inv_X[1][3]*F_filter->N4_X;
    F_filter->a3_X = F_filter->M_inv_X[2][0]*F_filter->N1_X + F_filter->M_inv_X[2][1]*F_filter->N2_X
    +F_filter->M_inv_X[2][2]*F_filter->N3_X + F_filter->M_inv_X[2][3]*F_filter->N4_X;
    F_filter->a4_X = F_filter->M_inv_X[3][0]*F_filter->N1_X + F_filter->M_inv_X[3][1]*F_filter->N2_X
    +F_filter->M_inv_X[3][2]*F_filter->N3_X + F_filter->M_inv_X[3][3]*F_filter->N4_X;
    
    F_filter->a1_AE = F_filter->M_inv_AE[0][0]*F_filter->N1_AE + F_filter->M_inv_AE[0][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[0][2]*F_filter->N3_AE + F_filter->M_inv_AE[0][3]*F_filter->N4_AE;
    F_filter->a2_AE = F_filter->M_inv_AE[1][0]*F_filter->N1_AE + F_filter->M_inv_AE[1][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[1][2]*F_filter->N3_AE + F_filter->M_inv_AE[1][3]*F_filter->N4_AE;
    F_filter->a3_AE = F_filter->M_inv_AE[2][0]*F_filter->N1_AE + F_filter->M_inv_AE[2][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[2][2]*F_filter->N3_AE + F_filter->M_inv_AE[2][3]*F_filter->N4_AE;
    F_filter->a4_AE = F_filter->M_inv_AE[3][0]*F_filter->N1_AE + F_filter->M_inv_AE[3][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[3][2]*F_filter->N3_AE + F_filter->M_inv_AE[3][3]*F_filter->N4_AE;
}

void get_Fstat_logL(struct Orbit *orbit, struct Data *data, double f0, double fdot, double theta, double phi, double *logL_X, double *logL_AE, double *Fparams)
{
    long M_filter, N_filter;
    long q;
    
    q = (long)(f0*data->T); 	// carrier frequency bin
    
    struct Filter *F_filter = malloc(sizeof(struct Filter));
    
    
    F_filter->f0     = f0;
    F_filter->fdot   = fdot;
    F_filter->fddot  = 0.;    //11.0/3.0*fdot*fdot/f0 
    F_filter->q      = q;
    F_filter->theta  = theta;
    F_filter->phi    = phi;
    

    int BW = galactic_binary_bandwidth(orbit->L, orbit->fstar, f0, fdot, cos(theta), 1.e-22, data->T, data->N);
    M_filter = BW;
    N_filter = BW;

    F_filter->M_filter = M_filter;
    F_filter->N_filter = N_filter;
    init_A_filters(orbit, data, F_filter);
    
    /////////
    //
    // Calculate N^{i} (Cornish & Crowder '05)
    // N^{i} = (s|A^{i})
    //
    /////////	
    get_N(data, F_filter);
    
    /////////
    //
    // Calculate M^{ij} inverse (Cornish & Crowder '05)
    // M^{ij} = (A^{i}|A^{j})
    // It is symmetric
    //
    /////////
    init_M_matrix(F_filter, data);	
    
    calc_Fstat_logL(F_filter);
    
    *logL_X  = F_filter->Fstat_X;
    *logL_AE = F_filter->Fstat_AE;
    
    calc_a_i(F_filter);           // calculate the a_{i}'s associated with filters to get F-stat ML params 
    get_F_params(F_filter);  // get the F-stat ML parameters	
    
    Fparams[0] = F_filter->A_AE_Fstat;
    Fparams[1] = F_filter->iota_AE_Fstat;
    Fparams[2] = F_filter->psi_AE_Fstat;
    Fparams[3] = F_filter->phase_AE_Fstat;
    
    free_Filter(F_filter);
}

void get_Fstat_xmax(struct Orbit *orbit, struct Data *data, double *x, double *xmax)
{
    struct Filter *F_filter = malloc(sizeof(struct Filter));

    double f0, theta, phi, fdot;
    
    {
        struct Source temp;
        map_array_to_params(&temp, x, data->T);

        f0 = temp.f0;
        theta = acos(temp.costheta);
        phi = temp.phi;
        fdot = temp.dfdt;
    }

    F_filter->f0     = f0;
    F_filter->fdot   = fdot;
    F_filter->fddot  = 0.;    //11.0/3.0*fdot*fdot/f0
    F_filter->q      = (long)x[F0];
    F_filter->theta  = theta;
    F_filter->phi    = phi;
    
 
    long M_filter, N_filter;
    
    int BW = galactic_binary_bandwidth(orbit->L, orbit->fstar, f0, fdot, cos(theta), 1.e-22, data->T, data->N);
    M_filter = BW;
    N_filter = BW;

    F_filter->M_filter = M_filter;
    F_filter->N_filter = N_filter;

    
    init_A_filters(orbit, data, F_filter);
    
    /////////
    //
    // Calculate N^{i} (Cornish & Crowder '05)
    // N^{i} = (s|A^{i})
    //
    /////////
    get_N(data, F_filter);
    
    /////////
    //
    // Calculate M^{ij} inverse (Cornish & Crowder '05)
    // M^{ij} = (A^{i}|A^{j})
    // It is symmetric
    //
    /////////
    init_M_matrix(F_filter, data);
    
    calc_a_i(F_filter);           // calculate the a_{i}'s associated with filters to get F-stat ML params
    get_F_params(F_filter);  // get the F-stat ML parameters
    
    xmax[AMP] = log(F_filter->A_AE_Fstat);
    xmax[COSI] = cos(F_filter->iota_AE_Fstat);
    xmax[PSI] = F_filter->psi_AE_Fstat;
    xmax[PHI0] = F_filter->phase_AE_Fstat;
    
    free_Filter(F_filter);
}



