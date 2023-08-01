#include <stdio.h>  // for printf errors
#include <stdlib.h> // for exit(1)

#include "gb_params.h"
// The overall goals of this subsytem:
//
// 1) Make it easy for leaf calculation code to tell 
//    whether a parameter is in use, and reference parameters by name.
// 2) Decouple the choice of parameterization from the memory layout of double * params
// 3) Permit runtime flag based decisions about the number, type, and order of 
//    parameters without unnecessary changes to calculation code.
// 4) Store parameters contiguously in memory so for loops over all parameters are performant 
//    and vectorizable by the compiler


// These are the singleton state for setup_parameters
// p_Ext2Int being initialized to p_Ext2Int[i] = i supports the use of
// label defines in the initial call to setup_parameters.
//
// Note: When adding a parameter it is important that these two
// lines remained updated with respect to the total number of parameters
// allowed to be referenced in code.
int p_Ext2Int [] = {0,1,2,3,4,5,6,7,8,9,10,11};

// This is the size of p_Ext2Int. It should not need to be referenced
// outside of this file
#define PARAM_MAX (12)

// This is the number of parameters in use at runtime.
// zero is reserved as a flag value meaning that setup_parameters
// has not yet run.
int NUM_PARAMS = 0;

// TODO: Add a lookup table for label macros to string names. This will 
//    enable things like comments at the top of data files with 
//    names of parameters, and more ergonomic error messages/debug 
//    printfs within calculation code.
void setup_parameters(int * parameterList, int NP) {
	if(NUM_PARAMS!=0) {
		fprintf(stderr, "Error: Multiple calls to setup_parameters. Only the first call is valid. Exiting now.\n");
		exit(1);
	}

	if(NP>PARAM_MAX) {
		fprintf(stderr, "Error: setup_parameters called with %d parameters,\n", NP); 
		fprintf(stderr, "but only %d parameters are supported by this build.\n", PARAM_MAX - 1);
		fprintf(stderr, "See gb_params.h/c for definitions. Exiting now.\n");
		exit(1);
	}

	if(NP==0) {
		fprintf(stderr, "Error: setup_parameters called with zero parameters. Exiting now.\n");
		exit(1);
	}
	
	// Populate NOT_IN_USE values so is_param functions properly
	for(int i = 0; i < PARAM_MAX; i++){
		p_Ext2Int[i] = PARAM_NOT_IN_USE;
	}

	// Populate actual parameter labels so parameter label macros function.
	for(int j = 0; j < NP; j++) {
		p_Ext2Int[parameterList[j]] = j;
	}

	// Set the global NUM_PARAMS so it can be used for iteration
	NUM_PARAMS = NP;
}

