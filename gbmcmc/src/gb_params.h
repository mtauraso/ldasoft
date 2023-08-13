
#ifndef _GB_PARAMS_H
#define _GB_PARAMS_H

#include <limits.h>   // for INT_MAX

// This is an internal index translation table which converts
// From the parameter label macros below (F0, COSTHETA, etc)
// to the actual index that params is laid out with.
//
// Code outside of params.c should never alter this directly
// only use it through the below macros
extern int p_Ext2Int [];

// After setup_parameters is run this contains the number of
// parameters for this run, and can be used for loops over
// all parameters.
//
// Example: Iterate over all parameters except amplitude
//
// 	for(int i=0; i < NUM_PARAMS; i++) {
// 		if (i == AMP) continue;
// 		params[i] = doMath(params[i]);
// 	}
//
// Example: Allocate a paramaters array
//
// 	double * params = calloc(NUM_PARAMS, 0)
//
extern int NUM_PARAMS;

// These are macros for use with setup_parameters and with 
// indexing the params array at runtime. Regardless of the 
// layout set up by setup_parameters params[F0] will always
// give you the principle frequency parameter, and so on.
// 
// Note: when adding a parameter params.c must also be updated
//       because it contains the length of this translation
//       array and initial values.
//
// TODO: These are accessing individual bits if there are cache/dereference
//       performance issues in leaf code, these could be converted to a single
//       int p_Ext2Int and accessed as a bitfield. 
#define F0             (p_Ext2Int[0])
#define F0_DESC        "Principle frequency unitless ratio [Hz*(data duration)^-1]"
#define COSTHETA       (p_Ext2Int[1])
#define COSTHETA_DESC  "Sky location Cos of latitude angle in elliptic coordinates dimensionless [-1,1]"
#define PHI            (p_Ext2Int[2])
#define PHI_DESC       "Sky location longitude angle in elliptic coordinates radians [0, 2pi]"
#define AMP            (p_Ext2Int[3])
#define AMP_DESC       "Gravitational wave amplitude log(strain)" // TODO Are units right?
#define COSI           (p_Ext2Int[4])
#define COSI_DESC      "Inclination Cos of inclination angle to ecliptic [-1,1]" // TODO verify, put range and units
#define PSI            (p_Ext2Int[5])
#define PSI_DESC       "Polarization angle. Radians [-pi,pi]" // TODO units and range
#define PHI0           (p_Ext2Int[6])
#define PHI0_DESC      "Phase of the binary at observation start radians [0, 2pi]" // Todo units and range
#define DFDT           (p_Ext2Int[7])
#define DFDT_DESC      "First derivative of principle frequency of gravitational wave. unitless ratio [Hz * (data duration)^-2]"
#define D2FDT2         (p_Ext2Int[8])
#define D2FDT2_DESC    "Second derivative of principle frequency of gravitational wave. unitless ratio [Hz * (data duration)^-3]"
// Above this line are legacy parameters. If you see a line or comment that says 
// params[x] where x is some constant, look at the indicies above for its name
// e.g. params[3] in the old system was AMP.
// Below this line there should no potential sensitivity of legacy code to the index numbers used.
#define DIST           (p_Ext2Int[9])
#define DIST_DESC      "Luminosity distance of source. Measured in kpc"
#define MC             (p_Ext2Int[10])
#define MC_DESC        "Chirp mass of binary pair in M_sun units"
#define DFDTASTRO      (p_Ext2Int[11])
#define DFDTASTRO_DESC "First derivative of principle frequency of binary system not attributable to gravitational waves. unitless ratio [Hz * (data duration)^-2]"

#define PARAM_NOT_IN_USE     (INT_MAX)


// This function defines a layout for the parameters in memory
// parameterList is a list of paremeter label defines, and NP 
// is the length of the parameter labels passed
// 
// IMPORTANT: This must be run before any thing else in params.h can be used. 
//            This can only run only once.
//
// Example usage:
// With a static set of parameters:
//
// 	int parameterList[] = {F0, COSTHETA, PI}
// 	setup_parameters(parameterList, 3);
//
// With a dynamic set of parameters dependent on a flag.
// 
// 	int NP = flags->skyLocation ? 3 : 1;
//      int * parameterList = malloc(NP * sizeof(int));
//      parameterList[0] = F0;
//      if(flags->skyLocation) {
//      	parameterList[1] = COSTHETA;
//      	parameterList[2] = PHI;
//      }
//      setup_parameters(parameterList, NP);
//      free(parameterList);
//
void setup_parameters(int * parameterList, int NP);

// This is only valid to call after setup_parameters, and tells
// you whether a particular paremeter is part of the layout
// for this run. 
//
// This is most useful when you need a parameter that could
// be a base parameter or a derived parameter depending on 
// runtime flags.
//
// int amp;
// if (is_param(AMP)) {
// 	amp = exp(params[AMP]);
// } else {
// 	amp = deriveAmplitude(params);
// }
#define is_param(label) (label != PARAM_NOT_IN_USE)

#endif // _GB_PARAMS_H

