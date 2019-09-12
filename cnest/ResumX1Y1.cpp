#include "ResumEFT.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#define LIMIT 10000
#define EPSABS 0
#define EPSREL 1e-4
#define MIN 0.0001
#define MAX 1


double X1 (const double & q, const InterpFunc & p) {
	return 2.6 * gsl_spline_eval (p.interp, q, p.accel) ;
}

double Y1 (const double & q, const InterpFunc & p) {
	return gsl_spline_eval (p.interp, q, p.accel) ;
}


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////

struct ParamsIntegr_ResumX1 {
	double q ;
	ParamsP11 InterpP11 ;
} ;

double ResumX1_Integrand_GSL (double x, void * params) {
	ParamsIntegr_ResumX1 p = *(ParamsIntegr_ResumX1 *) params ;

	return exp(-pow(x/CutResum,2)) * P11(x, p.InterpP11) * ( 2./3. -2.*gsl_sf_bessel_j1(x*p.q)/(x*p.q) ) ;
}

double ResumX1 (const double & q, const ParamsP11 & InterpP11) {

	ParamsIntegr_ResumX1 p = { q, InterpP11 } ;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (LIMIT);
	gsl_function F ;
  	F.function = &ResumX1_Integrand_GSL ;
  	F.params = &p ;

  	double res, err ;
	gsl_integration_qag (&F, MIN, MAX, EPSABS, EPSREL, LIMIT, GSL_INTEG_GAUSS41, w, &res, &err) ; 

	gsl_integration_workspace_free (w) ;

	return 0.5/pow(M_PI,2) * res ;
}



double ResumY1_Integrand_GSL (double x, void * params) {
	ParamsIntegr_ResumX1 p = *(ParamsIntegr_ResumX1 *) params ;

	return exp(-pow(x/CutResum,2)) * P11(x, p.InterpP11) * ( -2.*gsl_sf_bessel_j0(x*p.q) + 6.*gsl_sf_bessel_j1(x*p.q)/(x*p.q) ) ;
}

double ResumY1 (const double & q, const ParamsP11 & InterpP11) {

	ParamsIntegr_ResumX1 p = { q, InterpP11 } ;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (LIMIT);
	gsl_function F ;
  	F.function = &ResumY1_Integrand_GSL ;
  	F.params = &p ;

  	double res, err ;
	gsl_integration_qag (&F, MIN, MAX, EPSABS, EPSREL, LIMIT, GSL_INTEG_GAUSS41, w, &res, &err) ; 

	gsl_integration_workspace_free (w) ;

	return 0.5/pow(M_PI,2) * res ;
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////

void GetX1Y1 (const ParamsP11 & InterpP11, InterpFunc & InterpX1, InterpFunc & InterpY1, double & X1Infinity) {
	
/*	double X1Table[Nkp], Y1Table[Nkp] ;

	for (unsigned int m = 0 ; m < Nkp ; m++) {
		X1Table[m] = ResumX1(kplist[m],InterpP11) ;
		Y1Table[m] = ResumY1(kplist[m],InterpP11) ;
	}

	InterpX1.accel = gsl_interp_accel_alloc () ;
	InterpX1.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
	gsl_spline_init (InterpX1.interp, kplist, X1Table, Nkp) ;

	InterpY1.accel = gsl_interp_accel_alloc () ;
	InterpY1.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
	gsl_spline_init (InterpY1.interp, kplist, Y1Table, Nkp) ;

	X1Infinity = ResumX1(qInfinity,InterpP11) ;*/
}


