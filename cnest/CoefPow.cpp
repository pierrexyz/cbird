#include "CBiRd.h"
#include <fftw3.h>

void CoefPow (dcomplex Coef[], dcomplex Pow[], const ParamsP11 & paramsP11, const double & bias, const size_t & Nmax, const double & kmin, const double & kmax) {
	
	double dk = log(kmax/kmin) / ((double)Nmax-1.) ;
	
	double k[Nmax], Pn[Nmax] ;
	for (unsigned int i = 0 ; i < Nmax ; i++) {
		k[i] = kmin * exp(i*dk) ;
		Pn[i] = P11(k[i],paramsP11) * exp(-bias*i*dk) ;
	}

	for (unsigned int m = 0 ; m < Nmax+1 ; m++) Pow[m] = { bias, 2.*M_PI / (Nmax*dk) * ((double)m-(double)Nmax/2.) } ;

	dcomplex tmp[Nmax/2+1] ;
	fftw_plan forward_plan = fftw_plan_dft_r2c_1d(Nmax, Pn, (fftw_complex *) tmp, FFTW_ESTIMATE) ;
	fftw_execute(forward_plan) ;
	fftw_destroy_plan(forward_plan) ; 

	for (unsigned int m = 0 ; m < Nmax+1 ; m++) {
		if (m < Nmax/2) Coef[m] = conj(tmp[Nmax/2-m]) * pow(kmin, -Pow[m]) / ((double) Nmax) ;
		else Coef[m] = tmp[m-Nmax/2] * pow(kmin, -Pow[m]) / ((double) Nmax) ;
	}

	Coef[0] /= 2. ; Coef[Nmax] /= 2. ;
} 

/* std=c++14
void CoefPow (dcomplex Coef[], dcomplex Pow[], const ParamsP11 & paramsP11, const double & bias, const size_t & Nmax, const double & kmin, const double & kmax) {
	
	double dk = log(kmax/kmin) / ((double)Nmax-1.) ;
	
	double k[Nmax], Pn[Nmax] ;
	for (unsigned int i = 0 ; i < Nmax ; i++) {
		k[i] = kmin * exp(i*dk) ;
		Pn[i] = P11(k[i],paramsP11) * exp(-bias*i*dk) ;
	}

	for (unsigned int m = 0 ; m < Nmax+1 ; m++) Pow[m] = bias + 1i * 2.*M_PI / (Nmax*dk) * ((double)m-(double)Nmax/2.) ;

	dcomplex tmp[Nmax/2+1] ;
	fftw_plan forward_plan = fftw_plan_dft_r2c_1d(Nmax, Pn, (fftw_complex *) tmp, FFTW_ESTIMATE) ;
	fftw_execute(forward_plan) ;
	fftw_destroy_plan(forward_plan) ; 

	for (unsigned int m = 0 ; m < Nmax+1 ; m++) {
		if (m < Nmax/2) Coef[m] = conj(tmp[Nmax/2-m]) * pow(kmin, -Pow[m]) / ((double) Nmax) ;
		else Coef[m] = tmp[m-Nmax/2] * pow(kmin, -Pow[m]) / ((double) Nmax) ;
	}

	Coef[0] /= 2. ; Coef[Nmax] /= 2. ;
} 
*/
