#include "CBiRd.h"
#include "ResumEFT.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#define EPSABS 0
#define EPSREL 1e-5

double Mtest (const double & q, const InterpFunc & p) {
	return gsl_spline_eval (p.interp, q, p.accel) ;
}

double ResumM_Integrand_GSL (double q, void * params) {
	ParamsResumIntegr p = *(ParamsResumIntegr *) params ;

	return 0.5*pow(q/M_PI,2) *Mtest(q,p.InterpM);
}

// Resummation of 1 * M matrices: should give 1 for all diagonal components M_00, M_22, etc. and 0 for all off diagonal, M_02, M_20, etc.
void ResumSanityCheck(const string & PathToOutput, const ParamsP11 & paramsP11, StoreM * TableM) {

	// Create M matrices and export them
	ResumM (PathToOutput, paramsP11, true, TableM) ;

	ParamsResumIntegr integr ;
	double X1qInfinity = ResumX1(qInfinity,paramsP11) ;

	// gsl integration stuff
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000) ;
	gsl_function F ;
  	F.function = &ResumM_Integrand_GSL ;

	for (unsigned order = 0 ; order < 2 ; order++) {

		unsigned int Morder, Nlp ;

		// Linear power spectra
		if (order == 0) {
			Morder = 1 ;
			Nlp = 3 ;
		}

		// 1-loop power spectra
		if (order == 1) {
			Morder = 0 ;
			Nlp = 3 ; // Nlp = 5 ; 
		}

		double integrM[Nlout][Nlp][Nout] ;

		for (unsigned int i = 0 ; i < Nlout ; i++) { // loop over l = 0,2,4
			unsigned l = 2*i ;
			for (unsigned int m = 0 ; m < Nout ; m++) { // loop over k
				double k = kout[m] ;
				for (unsigned int j = 0 ; j < Nlp ; j++) { // loop over lp
					unsigned lp = 2*j ;

					integrM[i][j][m] = 0. ;

					integr.InterpM.accel = gsl_interp_accel_alloc () ;
					integr.InterpM.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
					gsl_spline_init (integr.InterpM.interp, (*TableM)[0][Morder][i][j][m], (*TableM)[1][Morder][i][j][m], Nkp) ;

					F.params = &integr ;

					double res = 0, err = 0 ;

					gsl_set_error_handler_off () ;
					gsl_integration_qag (&F, CutIRresum, CutUVresum, EPSABS, EPSREL, 100000, GSL_INTEG_GAUSS61, w, &res, &err) ;
					
					double Qinf = 0. ;
					Qinf = QInfinity(Morder,k,l,lp,X1qInfinity, paramsP11.f) ;

					integrM[i][j][m] = res + Qinf ;

					UnloadInterp (integr.InterpM) ;
				}
			}
		}

		// Output in files:
		int pre = 16 ; // white spacing 

		for (unsigned int i = 0 ; i < Nlout ; i++) { // loop over l = 0,2,4
			unsigned l = 2*i ;
			for (unsigned int j = 0 ; j < Nlp ; j++) { // loop over lp
				unsigned lp = 2*j ;
				ostringstream filename ;
				filename << PathToOutput << "/resum_data/integrM" << Morder << "_" << l << lp << ".dat" ;
				ofstream write (filename.str(), ios::out | ios::trunc) ;
		    	for (unsigned int m = 0 ; m < Nout ; m++) write << setw(pre) << kout[m] << " " << setw(pre-1) << integrM[i][j][m] << endl ;
				write.close() ;
			}
		}
	}

	gsl_integration_workspace_free (w) ;
}