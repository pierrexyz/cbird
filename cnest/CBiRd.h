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

#include "klist.h"

using namespace std ;

typedef complex<double> dcomplex ;

// Stuffs for FFTLog decomposition of the linear power spectrum
const size_t NFFT = 256 ;
const double kminFFT = 1.5e-5 ;
const double kmaxFFT = 20. ;
const double biasFFT = -1.6 ;
const double k0 = 1e-3 ;

const size_t N22 = 20 ;
const size_t N13 = 8 ;

// IR-resummation cutoff
const double CutResum = 0.066;

// Math stuff
double Heaviside (const double & a) ;
const double Pi = M_PI ;
const double E = exp(1.) ;

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

typedef double redshift ;

/** Cosmology **/
const size_t Nc = 5 ; // Number of cosmological parameters 
typedef double ParametersCosmology[Nc] ; // cosmological parameters
const string ParametersCosmologyNames[Nc] = { "A_s", "n_s", "h", "omega_b", "omega_cdm" } ; // A_s refers to ln(10^10 A_s)

/* Reference cosmology: Planck2015 */
const ParametersCosmology Reference = { 3.094, 0.9645, 0.6727, 0.04917, 0.2647 } ;

// Linear Growth rate f: GrowthFunction.cpp
double GrowthFactor (double & Omega_m, const double & a) ;
double LinearGrowthRate (const ParametersCosmology & p, const redshift & z) ;

// LoadConfigFile.cpp
typedef bool YesNo ;
void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, 
	string & PathToFolder, string & PathToLinearPowerSpectrum, 
	YesNo & ComputePowerSpectrum, YesNo & ResumPowerSpectrum, YesNo & ComputeBispectrum,
	string & PathToTriangles, double & aperp, double & apar) ;

// LinearPowerSpectrum.cpp
void LoadP11 (const string & LinearPowerSpectrumData, const ParametersCosmology & cosmo, const redshift & z, ParamsP11 & p) ;
double P11 (const double & q, const ParamsP11 & params) ;
void UnloadP11 (const ParamsP11 & params) ;

// EisensteinHuPowerSpectrum.cpp
void LoadEinseinsteinHu(const ParametersCosmology & cosmo, const redshift & z0, ParamsP11 & params, const size_t & Nmax = NFFT, const double & kmin = kminFFT, const double & kmax = kmaxFFT) ;

// Multipole moments
const size_t Nl = 5 ;
typedef double MultipoleMoments[Nl] ; // l = 0, 2, 4, 6, 8
// We call mi with i: power of mu
const MultipoleMoments m0 = { 1., 0., 0., 0., 0. } ;
const MultipoleMoments m2 = { 1./3., 2./3., 0., 0., 0. } ;
const MultipoleMoments m4 = { 1./5., 4./7., 8./35., 0., 0. } ;
const MultipoleMoments m6 = { 1./7., 10./21., 24./77., 16./231., 0. } ;
const MultipoleMoments m8 = { 1./9., 40./99., 48./148., 64./495., 128./6435. } ;

// For each loop terms we associate a pointer to multipole moments.
const MultipoleMoments * const Multipole22[N22] = { &m2, &m2, &m2, &m4, &m0, &m0, &m0, &m0, &m0, &m0, &m2, &m4, &m6, &m2, &m4, &m2, &m4, &m4, &m6, &m8 } ;
const MultipoleMoments * const Multipole13[N13] = { &m2, &m0, &m0, &m2, &m4, &m2, &m4, &m6 } ;


// ComputePowerSpectra.cpp

const size_t N0 = 3 ; 	// 3 linear terms
//const size_t N1 = 21 ; // 12 1-Loop PowerSpectra + 9 CounterTerms
const size_t N1 = 24 ;//18 ; // 12 1-Loop PowerSpectra + 6 CounterTerms ; no more support for stochastic term
typedef double PowerSpectraNoResum[Nl][N1][Nk] ;

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps) ;
void ComputePowerSpectra1LoopNoResum (const ParamsP11 & params, const double & nbar, const double & km, const double & knl, PowerSpectraNoResum * Ps) ;

// ComputeLoopIntegrals.cpp
void ComputeLoopIntegrals (const ParamsP11 & params, PowerSpectraNoResum * Ps) ;

// ComputeCounterTerms.cpp
void ComputeCounterTerms (const ParamsP11 & params, const double & Nbar, const double & kM, const double & kNL, PowerSpectraNoResum * Ps) ;

// CoefPow.cpp
void CoefPow (dcomplex Coef[], dcomplex Pow[], const ParamsP11 & p, const double & bias = biasFFT, const size_t & Nmax = NFFT, const double & kmin = kminFFT, const double & kmax = kmaxFFT) ;

// M22.cpp
dcomplex I(const dcomplex & n1, const dcomplex & n2) ;
dcomplex M22 (const unsigned int & j, const dcomplex & n1, const dcomplex & n2, const dcomplex & In1n2, const double & f1) ;

// M13.cpp
dcomplex M13 (const unsigned int & j, const dcomplex & n1, const double & f1) ;

// ExportPowerSpectraNoResum.cpp
void ExportPowerSpectraNoResum (const string & PathToFolder, const unsigned & order, PowerSpectraNoResum * Ps) ;



///////////// TREE-LEVEL BISPECTRUM MONOPOLE WITH AP EFFECT //////////////////////
struct ParamsIntegrBispectrumAP {
	double k1,k2,k3 ;
	ParamsP11 p11 ;
	double nbar ;
	double aperp, apar ;
	int id ;
} ;

// IntegrandBispectrumAP.cpp
double IntegrandBispectrumAP (const int & id, const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, const ParamsP11 & p, const double & nbar, const double & aperp, const double & apar) ;

void ScoccimarroTransform (const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) ;

void ScoccimarroTransformWithAP (const double & aperp, const double & apar, 
	const double & k1, const double & k2, const double & k3, const double & mu1, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) ;

// ComputeBispectrum.cpp
void ComputeBispectrumMonopole (const string & PathToFolder, const string & PathToTriangles, const ParamsP11 & params, const double & nbar, const double & aperp, const double & apar) ;
static int IntegrandBispectrumAP_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) ;



#endif
