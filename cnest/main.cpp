#include "CBiRd.h"

#include <ctime>
#include <ratio>
#include <chrono>
using namespace std::chrono;

static Correlator Cf ;
static Correlator Ps ;
static Correlator Cfsmooth ;

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		cerr << "Error: no configuration file specified." << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	else {
		// Default values
		double nbar = 3e-4 , knl = 1., km = 1. ;
		redshift z0 = 0.61 ;
		ParametersCosmology cosmo ; for (unsigned int i = 0 ; i < Nc ; i++) cosmo[i] = Reference[i] ;
		YesNo ExportPowerSpectraNoResum = 1, ExportCorrFuncNoResum = 1, ResumPowerSpectrum = 1 ;
		string PathToOutput = "./", PathToLinearPowerSpectrum ;
		double aperp = 1., apar = 1. ;
		unsigned int Nlout = 2 ;

		LoadConfigFile (argv[1], nbar, km, knl, z0, cosmo, 
			PathToOutput, PathToLinearPowerSpectrum, 
			ExportPowerSpectraNoResum, ExportCorrFuncNoResum, ResumPowerSpectrum, 
			aperp, apar, Nlout) ;

		knl = km ;
		///////////////////////////
		//
		///////////////////////////

		ParamsP11 paramsP11 ;
		LoadP11 (PathToLinearPowerSpectrum, cosmo, z0, paramsP11) ;

		ParamsP11 paramsP11smooth ; 
		// if (cosmo[5] == 0) LoadEinseinsteinHu(cosmo, z0, paramsP11smooth) ;
		// else LoadEinseinsteinHuWithMassiveNu(cosmo, z0, paramsP11smooth) ;

		LoadEinseinsteinHu(cosmo, z0, paramsP11smooth) ;

		double s[Nx] ;
		size_t Nq = Nopti ;
		for (unsigned int i = 0 ; i < Nq ; i++) s[i] = qopti[i] ;

		// double s[Nx] ;
		// size_t Nq = Nopti+3 ;
		// for (unsigned int i = 0 ; i < 3 ; i++) s[i] = (1.+i) * 1e-2 ;
		// for (unsigned int i = 0 ; i < Nopti+3 ; i++) s[i+3] = qopti[i] ;

		double k[Nx] ;
		size_t Nk = Nout ;
		for (unsigned int i = 0 ; i < Nk ; i++) k[i] = kout[i] ;

		auto start = high_resolution_clock::now() ;

		ComputeCorrelator(paramsP11, paramsP11smooth, Nq, &s, &Cf, &Cfsmooth, Nk, &k, &Ps, Nlout) ;

		// for (unsigned int i = 0 ; i < 3 ; i++)
		// 	for (unsigned int l = 0 ; l < 3 ; l++)
		// 		for (unsigned int n = 0 ; n < N0+N1 ; n++) Cf[l][n][i] = 0. ;
		
		if (ExportPowerSpectraNoResum == true) ExportCorrelator(PathToOutput, "PowerSpectraNoResum", Nk, &k, &Ps, true, Nlout) ;

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop-start);
		cout << "CBiRd ran in " << duration.count() << " milliseconds." << endl ;

		ResumCorrelator (paramsP11, Nq, &s, &Cf, &Cfsmooth, &Ps, Nlout) ;

		start = high_resolution_clock::now();
		duration = duration_cast<milliseconds>(stop-start);
		cout << "CBiRd ran in " << duration.count() << " milliseconds." << endl ;

		ExportCorrelator(PathToOutput, "PowerSpectra", Nk, &k, &Ps, true, Nlout) ;
	}

	return 0 ;

}