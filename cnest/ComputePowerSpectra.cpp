#include "CBiRd.h"

void ComputePowerSpectra1LoopNoResum (const ParamsP11 & params, const double & nbar, const double & km, const double & knl, PowerSpectraNoResum * Ps) {
	ComputeLoopIntegrals (params, Ps) ;
	ComputeCounterTerms (params, nbar, km, knl, Ps) ;
}

///////

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps) {
	
	double f1 = params.f ;

	for (unsigned int i = 0 ; i < 3 ; i++) {
		for (unsigned int m = 0 ; m < Nk; m++) {
			double P11k = P11(klist[m],params) ;

			(*Ps)[i][0][m] = m4[i] *f1*f1 *P11k ; 	// 1
			(*Ps)[i][1][m] = m2[i] *2.*f1 *P11k ; 	// b1
			(*Ps)[i][2][m] = m0[i] *P11k ; 			// b1*b1
		}
	}
}
