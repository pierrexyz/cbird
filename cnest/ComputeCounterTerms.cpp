#include "CBiRd.h"

void ComputeCounterTerms (const ParamsP11 & params, const double & nbar, const double & km, const double & knl, PowerSpectraNoResum * Ps) {

	double km2 = km*km ;
	double knl2 = knl*knl ;

	double km4 = km*km*km*km ;
	double knl4 = knl*knl*knl*knl ;

	for (unsigned int i = 0 ; i < Nl ; i++) {
		for (unsigned int m = 0 ; m < Nk ; m++) {
			double P11k = P11(klist[m],params) ;
			double k2 = pow(klist[m],2) ;
			double k4 = pow(klist[m],4) ;
			(*Ps)[i][12][m] = m0[i] * 2./knl2 *k2 *P11k ;				// b1*b5
			(*Ps)[i][13][m] = m2[i] * 2./km2 *k2 *P11k ;				// b1*b6
			(*Ps)[i][14][m] = m4[i] * 2./km2 *k2 *P11k ; 				// b1*b7
			(*Ps)[i][15][m] = m2[i] * 2./knl2 *params.f *k2 *P11k ;		// b5
			(*Ps)[i][16][m] = m4[i] * 2./km2 *params.f *k2 *P11k ;		// b6
			(*Ps)[i][17][m] = m6[i] * 2./km2 *params.f *k2 *P11k ;		// b7

			(*Ps)[i][18][m] = m0[i] * 2./knl4 *k4 *P11k ;				// b1*b5
			(*Ps)[i][19][m] = m2[i] * 2./km4 *k4 *P11k ;				// b1*b6
			(*Ps)[i][20][m] = m4[i] * 2./km4 *k4 *P11k ; 				// b1*b7
			(*Ps)[i][21][m] = m2[i] * 2./knl4 *params.f *k4 *P11k ;		// b5
			(*Ps)[i][22][m] = m4[i] * 2./km4 *params.f *k4 *P11k ;		// b6
			(*Ps)[i][23][m] = m6[i] * 2./km4 *params.f *k4 *P11k ;		// b7
			//(*Ps)[i][18][m] = m0[i] * 1./nbar ;							// b8
			//(*Ps)[i][19][m] = m0[i] * 1./nbar/km2 *k2 ;					// b9
			//(*Ps)[i][20][m] = m2[i] * 1./nbar/km2 *params.f *k2 ; 		// b10
		}
	}
}

