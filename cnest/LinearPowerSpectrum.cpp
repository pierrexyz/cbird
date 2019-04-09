#include "CBiRd.h"

double P11 (const double & q, const ParamsP11 & p) {
	return gsl_spline_eval (p.interp, q, p.accel) ;
}

void LoadP11 (const string & LinearPowerSpectrumData, const ParametersCosmology & cosmo, const redshift & z0, ParamsP11 & params) {

	// Set f
	params.f = LinearGrowthRate (cosmo,z0) ;
	//cerr << "Linear growth rate:" << params.f << endl ;

	// Load Linear Power Spectrum P11 data 
	ifstream p11data(LinearPowerSpectrumData, ios::in) ;
	if (!p11data.is_open()) {
		cerr << "There was a problem opening the linear power spectrum data file:" << LinearPowerSpectrumData << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	size_t NPoints = 0 ;
	double ki, P11i ;
	
	// count points number
	while (p11data >> ki >> P11i) {
		if (NPoints == 0 && ki > kminFFT) {
			cerr << "Please choose a linear power spectrum data file with kmin < " << kminFFT << endl ;
			exit(EXIT_FAILURE) ;
		}
		NPoints++ ;
	}

	if (ki < kmaxFFT) { 
		cerr << "Please choose a linear power spectrum data file with kmax > " << kmaxFFT << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	// reset ifstream buffer at the beginning of the document
	p11data.clear() ;
	p11data.seekg(0, ios::beg) ;

	NPoints++ ;

	double kdata[NPoints], Plin[NPoints] ; 	 

	kdata[0] = 0. ;
	Plin[0] = 0. ;

	for (unsigned int i = 1 ; i < NPoints ; i++) {
		p11data >> kdata[i] >> Plin[i] ;
	}

	

	p11data.close() ;

	// Interpolate P11
	gsl_interp_accel * acc = gsl_interp_accel_alloc () ;
	gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, NPoints) ;
	gsl_spline_init (spline, kdata, Plin, NPoints) ;

	params.accel = acc ;
	params.interp = spline ;
}

void UnloadP11 (ParamsP11 & params) {
	gsl_spline_free (params.interp) ;
    gsl_interp_accel_free (params.accel) ;
}