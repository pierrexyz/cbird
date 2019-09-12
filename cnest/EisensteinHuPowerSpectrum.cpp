#include "CBiRd.h"

void LoadEinseinsteinHu(const ParametersCosmology & cosmo, const redshift & z0, ParamsP11 & params, const size_t & Nmax, const double & kmin, const double & kmax) {
  // Set f
  params.f = LinearGrowthRate (cosmo,z0) ;

  // Compute Einseinstein-Hu power spectrum
  double
    e         = exp(1),
    Theta_CMB = 2.7255/2.7,
    kstar     = 0.05,
    h         = cosmo[2] ,
    Omega_b   = cosmo[3]/h/h ,
    Omega_m   = (cosmo[3]+cosmo[4])/h/h,  
    ns        = cosmo[1],
    As        = exp(cosmo[0])*1e-10,
    D1        = GrowthFactor(Omega_m, 1./(1.+z0)) ;

  double dk = log(kmax/kmin) / ((double)Nmax-1.) ;
  
  double k[Nmax], Plin[Nmax] ;
  for (unsigned int i = 0 ; i < Nmax ; i++) {
    k[i] = kmin * exp(i*dk) ;
    // Transfer function
    double 
           alpha_Gamma  = 1 - 0.328*log(431*Omega_m*h*h) * Omega_b/Omega_m + 0.38*log(22.3*Omega_m*h*h) * pow(Omega_b/Omega_m,2),
           s            = 44.5 * log(9.83/(Omega_m*h*h)) / sqrt(1 + 10*pow(Omega_b*h*h, 0.75)) * h,
           Gamma_eff    = Omega_m*h * (alpha_Gamma + (1-alpha_Gamma)/(1 + pow(0.43*k[i]*s, 4))),
           q            = k[i] * pow(Theta_CMB, 2) / Gamma_eff,
           L0           = log(2*e + 1.8*q),
           C0           = 14.2 + 731/(1 + 62.5*q),
           T0           = L0/(L0 + C0*pow(q,2)) ;
    Plin[i] = 2*M_PI*M_PI * As * k[i] / pow(100/299792.458,4) * pow(k[i]/kstar, ns-1) * pow(T0*D1, 2) ;
  }

  // Interpolate P11
  gsl_interp_accel * acc = gsl_interp_accel_alloc () ;
  gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, Nmax) ;
  gsl_spline_init (spline, k, Plin, Nmax) ;

  params.accel = acc ;
  params.interp = spline ;
}