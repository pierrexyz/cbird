#include "ResumEFT.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define EPSABS 0
#define EPSREL 1e-5


//////////////////////////////////////////////////////////////////////

void ResumPowerSpectra (const string & PathToFolder, const ParamsP11 & params, PowerSpectraNoResum * Linear, PowerSpectraNoResum * Loop, const YesNo & ImportM, const YesNo & ExportM, StoreM * TableM) {

	ParamsResumIntegr integr ;
	double X1qInfinity = ResumX1(qInfinity,params) ;

	// gsl integration stuff
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000) ;
	gsl_function F ;
  	F.function = &Resum_Integrand_GSL ;


for (unsigned order = 0 ; order < 2 ; order++) {

	unsigned int Morder, Np, Nlp ;
	PowerSpectraNoResum * Ps ;

	// Linear power spectra
	if (order == 0) {
		Morder = 1 ;
		Np = N0 ;
		Nlp = 3 ;
		Ps = Linear ;
	}

	// 1-loop power spectra
	if (order == 1) {
		Morder = 0 ;
		Np = N1 ;
		Nlp = 3 ; // Nlp = 5 ; 
		Ps = Loop ;
	}
	
	double Pout[Nlout][Np][Nout] ; 

	for (unsigned int i = 0 ; i < Nlout ; i++) { // loop over l = 0,2,4

		unsigned l = 2*i ;

		for (unsigned int n = 0 ; n < Np ; n++) { // loop over Pn

			for (unsigned int m = 0 ; m < Nout ; m++) { // loop over k

				double k = kout[m] ;

				Pout[i][n][m] = 0. ;

				for (unsigned int j = 0 ; j < Nlp ; j++) { // sum over lp

					unsigned lp = 2*j ;

					// Interpolation of Pn
					integr.InterpP.accel = gsl_interp_accel_alloc () ;
					integr.InterpP.interp = gsl_spline_alloc (gsl_interp_cspline, Nk) ;
					gsl_spline_init (integr.InterpP.interp, klist, (*Ps)[j][n], Nk) ;

					if (ImportM == true) LoadM (PathToFolder,Morder,k,l,lp,integr.InterpM) ;

					else {
						integr.InterpM.accel = gsl_interp_accel_alloc () ;
						integr.InterpM.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
						gsl_spline_init (integr.InterpM.interp, (*TableM)[0][Morder][i][j][m], (*TableM)[1][Morder][i][j][m], Nkp) ;
					}

					F.params = &integr ;

					double res = 0, err = 0 ;

					gsl_set_error_handler_off () ;
					gsl_integration_qag (&F, CutIRresum, CutUVresum, EPSABS, EPSREL, 100000, GSL_INTEG_GAUSS61, w, &res, &err) ; 
					
					/* mid-point integration
					unsigned u = 0 ;
					while ((*TableM)[0][Morder][i][j][m][u] < 0.01) u++ ;

					while ((*TableM)[0][Morder][i][j][m][u] < 0.9) {
						u++ ;
						double k1 = (*TableM)[0][Morder][i][j][m][u] ;
						double k0 = (*TableM)[0][Morder][i][j][m][u-1] ;
						double M1 = (*TableM)[1][Morder][i][j][m][u] ;
						double M0 = (*TableM)[1][Morder][i][j][m][u-1] ;
						//cout <<  u << " " << k1 << endl  ;
						res +=  0.5 * ( 0.5*pow(k1/M_PI,2) *  P(k1,integr.InterpP) * M1 - 0.5*pow(k0/M_PI,2) *  P(k0,integr.InterpP) * M0 ) * (k1-k0) ;
						//cout << res << " " << k1-k0 << " " << M1 << endl ;
					}
					*/


					
					double Qinf = 0. ;
					if (ImportM == true) Qinf = LoadQInfinity(PathToFolder,Morder,k,l,lp) ;
					else Qinf = QInfinity(Morder,k,l,lp,X1qInfinity,params.f) ;

					Pout[i][n][m] += res + Qinf * P(k,integr.InterpP) ;
					

					//Pout[i][n][m] += res ; //+ Qinf * P(k,integr.InterpP) ;

					
					UnloadInterp (integr.InterpM) ;
					UnloadInterp (integr.InterpP) ;

				}
			}
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

//////////////////////////////////////////////////////////////////////

double Resum_Integrand_GSL (double q, void * params) {
	ParamsResumIntegr p = *(ParamsResumIntegr *) params ;

	return 0.5*pow(q/M_PI,2) * ( P(q,p.InterpP) *M(q,p.InterpM) ) ;
}

//////////////////////////////////////////////////////////////////////

double P (const double & q, const InterpFunc & p) {
	return gsl_spline_eval (p.interp, q, p.accel) ;
}

double M (const double & q, const InterpFunc & p) {
	return gsl_spline_eval (p.interp, q, p.accel) ;
}

//////////////////////////////////////////////////////////////////////
