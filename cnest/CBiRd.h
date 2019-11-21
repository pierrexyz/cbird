#ifndef CBIRD_H
#define CBIRD_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring> 
#include <sstream>
#include <fstream>
#include <complex>

#include <gsl/gsl_spline.h>

#include "Eigen/Dense"

#include "kout.h"

using namespace std ;

typedef complex<double> dcomplex ;

// min max IR-convolution
const double qMIN = 60. ;
const double qMAX = 170. ;

// Stuffs for FFTLog decomposition of the linear power spectrum
const size_t NFFT = 256 ;
const double kminFFT = 1.5e-5 ;
const double kmaxFFT = 20. ;
const double biasFFT = -1.6 ;
const double k0 = 1e-3 ;

const size_t N22 = 20 ;
const size_t N13 = 8 ;

// IR-resummation X1, Y1 cutoff
const double CutResum = 0.066 ;

/* Computes the Gamma function using the Lanczos approximation */
static dcomplex MyGamma(dcomplex z) {
    /* Lanczos coefficients for g = 7 */
    static double p[] = {
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
       -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
       -176.61502916214059906584551354,
        12.507343278686904814458936853,
       -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7
    };

    if(z.real() < 0.5)
        return M_PI / (sin(M_PI*z)*MyGamma(1. - z));

    z -= 1;
    dcomplex x = p[0];
    for(int n = 1; n < 9; n++)
        x += p[n] / (z + double(n));
    dcomplex t = z + 7.5;
    return sqrt(2*M_PI) * pow(t, z+0.5) * exp(-t) * x;
}

// Struct declarations
struct InterpFunc {
	gsl_interp_accel * accel ;
	gsl_spline * interp ;
} ;

struct ParamsP11 {
	double f ;
	gsl_interp_accel * accel ;
	gsl_spline * interp ;
} ;

struct ParamsResumIntegr {
    InterpFunc InterpCf0, InterpCf2, InterpCf4, InterpX1, InterpY1 ;
    int o, l, lp ;
    double k, f1 ;
} ;


/** Cosmology **/
const size_t Nc = 7 ; // Number of cosmological parameters 
typedef double ParametersCosmology[Nc] ; // cosmological parameters
const string ParametersCosmologyNames[Nc] = { "ln10^{10}A_s", "n_s", "h", "omega_b", "omega_cdm", "N_ncdm", "Sum_mnu" } ;

/* Reference cosmology: Planck2015 */
const ParametersCosmology Reference = { 3.094, 0.9645, 0.6727, 0.04917, 0.2647, 1, 0.06 } ;

// Linear Growth rate f: GrowthFunction.cpp
typedef double redshift ;
double GrowthFactor (double & Omega_m, const double & a) ;
double LinearGrowthRate (const ParametersCosmology & p, const redshift & z) ;

// LoadConfigFile.cpp
typedef bool YesNo ;
void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, 
	string & PathToOutput, string & PathToLinearPowerSpectrum, 
	YesNo & ExportPowerSpectrumNoResum, YesNo & ExportCorrFuncNoResum, YesNo & ResumPowerSpectrum,
	double & aperp, double & apar, unsigned int & Nlout) ;

// LinearPowerSpectrum.cpp
void LoadP11 (const string & LinearPowerSpectrumData, const ParametersCosmology & cosmo, const redshift & z, ParamsP11 & p) ;
double P11 (const double & q, const ParamsP11 & params) ;
void UnloadP11 (const ParamsP11 & params) ;

// EisensteinHuPowerSpectrum.cpp
void LoadEinseinsteinHu(const ParametersCosmology & cosmo, const redshift & z0, ParamsP11 & params, const size_t & Nmax = NFFT, const double & kmin = kminFFT, const double & kmax = kmaxFFT) ;
void LoadEinseinsteinHuWithMassiveNu(const ParametersCosmology & cosmo, const redshift & z0, ParamsP11 & params, const size_t & Nmax = NFFT, const double & kmin = kminFFT, const double & kmax = kmaxFFT) ;
int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm, int degen_hdm, float omega_lambda, float hubble, float redshift);
float TFmdm_onek_mpc(float kk) ;
float TFmdm_onek_hmpc(float kk) ;


// Multipole moments
const size_t Nl5 = 5 ;
const size_t Nl = 3 ;
typedef double MultipoleMoments[Nl5] ; // l = 0, 2, 4, 6, 8
// We call mi with i: power of mu
const MultipoleMoments m0 = { 1., 0., 0., 0., 0. } ;
const MultipoleMoments m2 = { 1./3., 2./3., 0., 0., 0. } ;
const MultipoleMoments m4 = { 1./5., 4./7., 8./35., 0., 0. } ;
const MultipoleMoments m6 = { 1./7., 10./21., 24./77., 16./231., 0. } ;
const MultipoleMoments m8 = { 1./9., 40./99., 48./148., 64./495., 128./6435. } ;

// For each loop terms we associate a pointer to multipole moments.
const MultipoleMoments * const Multipole22[N22] = { &m2, &m2, &m2, &m4, &m0, &m0, &m0, &m0, &m0, &m0, &m2, &m4, &m6, &m2, &m4, &m2, &m4, &m4, &m6, &m8 } ;
const MultipoleMoments * const Multipole13[N13] = { &m2, &m0, &m0, &m2, &m4, &m2, &m4, &m6 } ;

// CoefPow.cpp
void CoefPow (dcomplex Coef[], dcomplex Pow[], const ParamsP11 & p, const double & bias = biasFFT, const size_t & Nmax = NFFT, const double & kmin = kminFFT, const double & kmax = kmaxFFT) ;

// M22.cpp
dcomplex I(const dcomplex & n1, const dcomplex & n2) ;
dcomplex M22 (const unsigned int & j, const dcomplex & n1, const dcomplex & n2, const dcomplex & In1n2, const double & f1) ;

// M13.cpp
dcomplex M13 (const unsigned int & j, const dcomplex & n1, const double & f1) ;

// ComputeCorrelationFunction.cpp
const size_t N0 = 3 ;   // 3 linear terms
const size_t N1 = 18 ; // 12 1-Loop PowerSpectra + 6 CounterTerms ; no more support for stochastic term
const size_t Np = 21 ; 
const size_t Nx = 100 ; // max number of k or q 
typedef double Coordinates[Nx] ;
typedef double Correlator[Nl][N0+N1][Nx] ;
void ComputeCorrelator (const ParamsP11 & p, const size_t & Nq, Coordinates * q, Correlator * Cf , const size_t & Nk, Coordinates * k , Correlator * Ps, const unsigned int & Nlout = Nl) ;
dcomplex MPC (const unsigned int & l, const dcomplex & n1) ;

// ExportCorrelator.cpp
void ExportCorrelator (const string & Path, const char *, const size_t & Nmax, Coordinates * x, Correlator * O, const YesNo & ShotNoise = false, const unsigned int & Nlout = Nl) ;

// ResumQ.cpp
const double qInfinity = 10000. ;
const double alpha = 1./3. ;
const double beta = 1./3. ;
double ResumQ0 (const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & f1) ;
double ResumQ1(const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & f1) ;

// ResumX1Y1.cpp
double ResumX1 (const double & q, const ParamsP11 & InterpP11) ;
double ResumY1 (const double & q, const ParamsP11 & InterpP11) ;
void GetX1Y1 (const ParamsP11 & InterpP11, InterpFunc & InterpX1, InterpFunc & InterpY1, double & X1Infinity) ;
double X1 (const double & q, const InterpFunc & InterpX1) ;
double Y1 (const double & q, const InterpFunc & InterpX1) ;

// ResumCorrelator.cpp
double C (const double & q, const InterpFunc & InterpC) ;
double Q (const int & order, const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & ff) ;
void ResumCorrelator (const ParamsP11 & InterpP11, const size_t & Nq, Coordinates * q, Correlator * Cf, Correlator * Ps, const size_t & Nlout = Nl) ;
double Resum_Integrand_GSL (double q, void * params) ;
void UnloadInterp (InterpFunc & params) ;

static void ResizeArray(const size_t Na, double Array[], double SubArray[]) {
    for (unsigned int i = 0 ; i < Na ; i++) SubArray[i] = Array[i] ;
}

void ExtractBAO(const size_t & Nq, Coordinates * q, double Cf[], double BAO[]) ;

#endif
