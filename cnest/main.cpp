#include "CBiRd.h"

// #include <ctime>
// #include <ratio>
// #include <chrono>
// using namespace std::chrono;

static Correlator Cf ;
static Correlator Ps ;

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
		//LoadEinseinsteinHu(cosmo, z0, paramsP11) ;
		//LoadEinseinsteinHuWithMassiveNu(cosmo, z0, paramsP11) ;

		double s[Nx] ;
		size_t Nq = Nopti ;
		for (unsigned int i = 0 ; i < Nopti ; i++) s[i] = qopti[i] ;

		double k[Nx] ;
		size_t Nk = Nout ;
		for (unsigned int i = 0 ; i < Nk ; i++) k[i] = kout[i] ;

		// auto start = high_resolution_clock::now() ;

		ComputeCorrelator(paramsP11, Nq, &s, &Cf, Nk, &k, &Ps, Nlout) ;
		
		if (ExportPowerSpectraNoResum == true) ExportCorrelator(PathToOutput, "PowerSpectraNoResum", Nk, &k, &Ps, true, Nlout) ;

		ResumCorrelator (paramsP11, Nq, &s, &Cf, &Ps, Nlout) ;

		// start = high_resolution_clock::now();
		// duration = duration_cast<milliseconds>(stop-start);
		// cout << "CBiRd ran in " << duration.count() << " milliseconds." << endl ;

		ExportCorrelator(PathToOutput, "PowerSpectra", Nk, &k, &Ps, true, Nlout) ;
	}

	return 0 ;

}