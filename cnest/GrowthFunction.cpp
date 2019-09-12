#include "CBiRd.h"
#include <gsl/gsl_integration.h>

double HubbleParameter (const double & Omega_m, const double & a) {
	return sqrt( Omega_m * pow(a, -3) + 1. - Omega_m ) ;
}

double GrowthFactorIntegrand (double a, void * params) {	
	double Omega_m = *(double *) params ;
	return pow(a*HubbleParameter(Omega_m,a), -3) ;
}

double GrowthFactor (double & Omega_m, const double & a) {
	
	double H = HubbleParameter (Omega_m,a) ;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;
	
	double D, error ;
	
	gsl_function F ;
	F.function = &GrowthFactorIntegrand ;
	F.params = &Omega_m ;

	gsl_integration_qags (&F, 0, a, 0, 1e-6, 1000, w, &D, &error) ; 
	
	gsl_integration_workspace_free (w) ;
	
	return H * D ; 
}

double LinearGrowthRate (const ParametersCosmology & p, const redshift & z) {

	double Omega_m = (p[3]+p[4])/pow(p[2],2) ;
	double a = 1./(1.+z) ;
	double H = HubbleParameter (Omega_m,a) ;

	return -1.5*Omega_m /pow(a,3)/pow(H,2) + 1./pow(a*H,2)/GrowthFactor(Omega_m,a) ;
}
