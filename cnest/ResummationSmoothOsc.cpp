#include "ResumEFT.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define EPSABS 0
#define EPSREL 1e-4

void ResumSmoothOscPowerSpectra (const string & PathToFolder, const ParamsP11 & params, PowerSpectraNoResum * Linear, PowerSpectraNoResum * Loop, PowerSpectraNoResum * LinearSmooth, PowerSpectraNoResum * LoopSmooth, StoreM * TableM) {
	ParamsResumIntegr integr ;
	double X1qInfinity = ResumX1(qInfinity,params) ;

	double ResumOsc = 0, err = 0 ;

	// gsl integration stuff
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000) ;
	gsl_set_error_handler_off () ;
  	gsl_function Osc ;
  	Osc.function = &Resum_Integrand_GSL ;

	for (unsigned order = 0 ; order < 2 ; order++) {

		unsigned int Morder, Np, Nlp ;
		PowerSpectraNoResum * Ps, * Psmooth ;

		// Linear power spectra
		if (order == 0) { Morder = 1 ; Np = N0 ; Nlp = 3 ; Ps = Linear ; Psmooth = LinearSmooth ; }
		// 1-loop power spectra
		if (order == 1) { Morder = 0 ; Np = N1 ; Nlp = 3 ; Ps = Loop ; Psmooth = LoopSmooth ; }

		double Pout[Nlout][Np][Nout] ;

		// Resummation of smooth part: Int M_l,lp(k,kp) * Pn_lp,smooth(kp) = Pn_l,smooth(k) + O(1e-3) for l = lp, and 0 otherwise.
		// since Pn_lp,smooth(kp) is kp-independent (slowly varying over the kernel where M_l,lp(k,kp) is nonzero), 
		// and using Int M_l,lp(k,kp) = 1. for l = lp, and 0 otherwise.

		// Resummation of oscillating part
		for (unsigned int i = 0 ; i < Nlout ; i++) { // loop over l = 0,2,4
			unsigned l = 2*i ;
			for (unsigned int n = 0 ; n < Np ; n++) { // loop over Pn

				// Interpolation of Pn_l,smooth
				InterpFunc InterpSmooth ;
				InterpSmooth.accel = gsl_interp_accel_alloc () ;
				InterpSmooth.interp = gsl_spline_alloc (gsl_interp_cspline, Nk) ;
				gsl_spline_init (InterpSmooth.interp, klist, (*Psmooth)[i][n], Nk) ;

				for (unsigned int m = 0 ; m < Nout ; m++) { // loop over k
					double k = kout[m] ;
					Pout[i][n][m] = P(k, InterpSmooth) ;
					for (unsigned int j = 0 ; j < Nlp ; j++) { // sum over lp
						unsigned lp = 2*j ;

						//if (l == lp) {
						// Interpolation of Pn_lp - Pn_lp,smooth
						double Posc[Nk] ;
						for (unsigned int ii = 0 ; ii < Nk ; ii++) Posc[ii] = (*Ps)[j][n][ii] - (*Psmooth)[j][n][ii] ;
						integr.InterpP.accel = gsl_interp_accel_alloc () ;
						integr.InterpP.interp = gsl_spline_alloc (gsl_interp_cspline, Nk) ;
						gsl_spline_init (integr.InterpP.interp, klist, Posc, Nk) ;

						// Interpolation of M_l,lp (k, kp)
						integr.InterpM.accel = gsl_interp_accel_alloc () ;
						integr.InterpM.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
						gsl_spline_init (integr.InterpM.interp, (*TableM)[0][Morder][i][j][m], (*TableM)[1][Morder][i][j][m], Nkp) ;

						Osc.params = &integr ;

						gsl_integration_qag (&Osc, CutIRresum, CutUVresum, EPSABS, EPSREL, 100000, GSL_INTEG_GAUSS61, w, &ResumOsc, &err) ; 
						double Qinf = QInfinity(Morder,k,l,lp,X1qInfinity,params.f) ;

						Pout[i][n][m] += ResumOsc + Qinf * P(k, integr.InterpP) ;

						UnloadInterp (integr.InterpP) ;
						UnloadInterp (integr.InterpM) ;
						
						//}
					}
				}

				UnloadInterp (InterpSmooth) ;
			}
		}

		// Output in files:
		int pre = 16 ; // white spacing 

		ostringstream filename ;
		if (order == 0) filename << PathToFolder << "/PowerSpectraLinear.dat" ;
		if (order == 1) filename << PathToFolder << "/PowerSpectra1loop.dat" ;

		ofstream write (filename.str(), ios::out | ios::trunc) ;

		if (order == 0) {
			write << "#" << setw(pre-1) << "k" << setw(pre) << "1" << setw(pre) << "b1" << setw(pre) << "b1*b1" << endl ;
		}

		if (order == 1) { write << "#" << setw(pre-1) << "k" <<  setw(pre) << "1" << setw(pre) << "b1"  << setw(pre) << "b2"  << setw(pre) << "b3"  << setw(pre) << "b4"  << setw(pre) << "b1*b1"  << setw(pre) << "b1*b2"  << setw(pre) << "b1*b3"  << setw(pre) << "b1*b4"  << setw(pre) << "b2*b2"  << setw(pre) << "b2*b4"  << setw(pre) << "b4*b4" 
			<< setw(pre) << "b1*b5"  << setw(pre) << "b1*b6"  << setw(pre) << "b1*b7"  << setw(pre) << "b5"  << setw(pre) << "b6"  << setw(pre) << "b7"  << endl ;
			//<< setw(pre) << "b8"  << setw(pre) << "b9"  << setw(pre) << "b10" << endl ; 
		}

		for (unsigned int i = 0 ; i < Nlout ; i++) {
	    	for (unsigned int m = 0 ; m < Nout ; m++) {
	    		write << setw(pre) << kout[m] ;
	    		for (unsigned int p = 0 ; p < Np ; p++) write << " " << setw(pre-1) << Pout[i][p][m] ;  // no precision lack in standard output accuracy 								
	    		write << endl ;
	    	}
		}
		write.close() ;
	}

	gsl_integration_workspace_free (w) ;
}

/* Not so smooth, so not so well resummed ...
#define EPSABS 0
#define EPSREL 1e-4

#define SIGMA 0.1823 // log(1.2)

double GaussianWeight(const double & k, const double & sigma) {
	return exp(-0.5*k*k/sigma/sigma) / sqrt(2.*M_PI) / sigma ;
}

double Smooth_Integrand_GSL (double logkp, void * params) {
	ParamsResumIntegr p = *(ParamsResumIntegr *) params ;

	return P(exp(logkp), p.InterpP) * GaussianWeight(logkp-log(p.k), SIGMA) ;
}

void ResumSmoothOscPowerSpectra (const string & PathToFolder, const ParamsP11 & params, PowerSpectraNoResum * Linear, PowerSpectraNoResum * Loop, StoreM * TableM) {
	/*
	ParamsResumIntegr integr ;
	double X1qInfinity = ResumX1(qInfinity,params) ;

	double ResumSmooth = 0, ResumOsc = 0, err = 0 ;

	// gsl integration stuff
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000) ;
	gsl_set_error_handler_off () ;
	gsl_function Smooth ;
  	Smooth.function = &Smooth_Integrand_GSL ;
  	gsl_function Osc ;
  	Osc.function = &Resum_Integrand_GSL ;

	for (unsigned order = 0 ; order < 2 ; order++) {

		unsigned int Morder, Np, Nlp ;
		PowerSpectraNoResum * Ps ;

		// Linear power spectra
		if (order == 0) { Morder = 1 ;Np = N0 ; Nlp = 3 ; Ps = Linear ; }
		// 1-loop power spectra
		if (order == 1) { Morder = 0 ; Np = N1 ; Nlp = 3 ; Ps = Loop ; }

		double Pout[Nlout][Np][Nout] ; 
		double Psmooth[Nlp][Np][Nk] ;

		// Resummation of smooth part: Int M_l,lp(k,kp) * Pn_lp,smooth(kp) = Pn_l,smooth(k) + O(1e-3) for l = lp, and 0 otherwise.
		// since Pn_lp,smooth(kp) is kp-independent (slowly varying over the kernel where M_l,lp(k,kp) is nonzero), 
		// and using Int M_l,lp(k,kp) = 1. for l = lp, and 0 otherwise.
		for (unsigned int j = 0 ; j < Nlp ; j++) { // loop over lp
			unsigned lp = 2*j ;
			for (unsigned int n = 0 ; n < Np ; n++) { // loop over Pn
				for (unsigned int m = 0 ; m < Nk ; m++) { // loop over kp: we need to evaluate the smooth part over kp since we will subtract it to the total to get the oscillating part
					integr.k = klist[m] ;

					// Interpolation of Pn_lp
					integr.InterpP.accel = gsl_interp_accel_alloc () ;
					integr.InterpP.interp = gsl_spline_alloc (gsl_interp_cspline, Nk) ;
					gsl_spline_init (integr.InterpP.interp, klist, (*Ps)[j][n], Nk) ;

					Smooth.params = &integr ;

					gsl_integration_qag (&Smooth, log(CutIRresum), log(CutUVresum), EPSABS, EPSREL, 100000, GSL_INTEG_GAUSS41, w, &ResumSmooth, &err) ; 
					Psmooth[j][n][m] = ResumSmooth ;

					UnloadInterp (integr.InterpP) ;
				}
			}
		}

		// Oscillating part
		for (unsigned int i = 0 ; i < Nlout ; i++) { // loop over l = 0,2,4
			unsigned l = 2*i ;
			for (unsigned int n = 0 ; n < Np ; n++) { // loop over Pn

				// Interpolation of Pn_l,smooth
				InterpFunc InterpSmooth ;
				InterpSmooth.accel = gsl_interp_accel_alloc () ;
				InterpSmooth.interp = gsl_spline_alloc (gsl_interp_cspline, Nk) ;
				gsl_spline_init (InterpSmooth.interp, klist, Psmooth[i][n], Nk) ;

				for (unsigned int m = 0 ; m < Nout ; m++) { // loop over k
					double k = kout[m] ;
					Pout[i][n][m] = P(k, InterpSmooth) ;
					for (unsigned int j = 0 ; j < Nlp ; j++) { // sum over lp
						unsigned lp = 2*j ;

						//if (l == lp) {
						
						integr.InterpP.accel = gsl_interp_accel_alloc () ;
						integr.InterpP.interp = gsl_spline_alloc (gsl_interp_cspline, Nk) ;
						
						// Interpolation of Pn_lp - Pn_lp,smooth
						double Posc[Nk] ;
						for (unsigned int ii = 0 ; ii < Nk ; ii++) Posc[ii] = (*Ps)[j][n][ii] - Psmooth[j][n][ii] ;
						gsl_spline_init (integr.InterpP.interp, klist, Posc, Nk) ;

						// Interpolation of M_l,lp (k, kp)
						integr.InterpM.accel = gsl_interp_accel_alloc () ;
						integr.InterpM.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
						gsl_spline_init (integr.InterpM.interp, (*TableM)[0][Morder][i][j][m], (*TableM)[1][Morder][i][j][m], Nkp) ;

						Osc.params = &integr ;

						gsl_integration_qag (&Osc, CutIRresum, CutUVresum, EPSABS, EPSREL, 100000, GSL_INTEG_GAUSS61, w, &ResumOsc, &err) ; 
						double Qinf = QInfinity(Morder,k,l,lp,X1qInfinity,params.f) ;

						Pout[i][n][m] += ResumOsc + Qinf * P(k, integr.InterpP) ;

						UnloadInterp (integr.InterpP) ;
						UnloadInterp (integr.InterpM) ;
						
						//}
					}
				}

				UnloadInterp (InterpSmooth) ;
			}
		}

		// Output in files:
		int pre = 16 ; // white spacing 

		ostringstream filename ;
		if (order == 0) filename << PathToFolder << "/PowerSpectraLinear.dat" ;
		if (order == 1) filename << PathToFolder << "/PowerSpectra1loop.dat" ;

		ofstream write (filename.str(), ios::out | ios::trunc) ;

		if (order == 0) {
			write << "#" << setw(pre-1) << "k" << setw(pre) << "1" << setw(pre) << "b1" << setw(pre) << "b1*b1" << endl ;
		}

		if (order == 1) { write << "#" << setw(pre-1) << "k" <<  setw(pre) << "1" << setw(pre) << "b1"  << setw(pre) << "b2"  << setw(pre) << "b3"  << setw(pre) << "b4"  << setw(pre) << "b1*b1"  << setw(pre) << "b1*b2"  << setw(pre) << "b1*b3"  << setw(pre) << "b1*b4"  << setw(pre) << "b2*b2"  << setw(pre) << "b2*b4"  << setw(pre) << "b4*b4" 
			<< setw(pre) << "b1*b5"  << setw(pre) << "b1*b6"  << setw(pre) << "b1*b7"  << setw(pre) << "b5"  << setw(pre) << "b6"  << setw(pre) << "b7"  << endl ;
			//<< setw(pre) << "b8"  << setw(pre) << "b9"  << setw(pre) << "b10" << endl ; 
		}

		for (unsigned int i = 0 ; i < Nlout ; i++) {
	    	for (unsigned int m = 0 ; m < Nout ; m++) {
	    		write << setw(pre) << kout[m] ;
	    		for (unsigned int p = 0 ; p < Np ; p++) write << " " << setw(pre-1) << Pout[i][p][m] ;  // no precision lack in standard output accuracy 
	    																																
	    		write << endl ;
	    	}
		}
		write.close() ;
	}

	gsl_integration_workspace_free (w) ;
}
*/

