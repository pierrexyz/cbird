#include "CBiRd.h"
#include "ResumEFT.h"

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
		YesNo ComputePowerSpectrum = 1, ResumPowerSpectrum = 1, ComputeBispectrum = 0 ;
		YesNo ExportM = 0, ImportM = 0 ;
		string PathToOutput = "./", PathToLinearPowerSpectrum, PathToTriangles ;
		double aperp = 1., apar = 1. ;

		LoadConfigFile (argv[1], nbar, km, knl, z0, cosmo, 
			PathToOutput, PathToLinearPowerSpectrum, 
			ComputePowerSpectrum, ResumPowerSpectrum, ComputeBispectrum,
			PathToTriangles, aperp, apar) ;

		knl = km ;
		///////////////////////////
		//
		///////////////////////////

		ParamsP11 paramsP11 ;
		LoadP11 (PathToLinearPowerSpectrum, cosmo, z0, paramsP11) ;

		if (ComputePowerSpectrum == true) {

			PowerSpectraNoResum Ps1Loop ;
			PowerSpectraNoResum PsLinear ;

			ComputePowerSpectraLinearNoResum  (paramsP11, &PsLinear) ;
			ComputePowerSpectra1LoopNoResum (paramsP11, nbar, km, knl, &Ps1Loop) ;
			
			if (ResumPowerSpectrum == true) {

				ParamsP11 paramsP11Smooth ;
				LoadEinseinsteinHu(cosmo, z0, paramsP11Smooth) ;

				PowerSpectraNoResum Ps1LoopSmooth ;
				PowerSpectraNoResum PsLinearSmooth ;

				ComputePowerSpectraLinearNoResum  (paramsP11Smooth, &PsLinearSmooth) ;
				ComputePowerSpectra1LoopNoResum (paramsP11Smooth, nbar, km, knl, &Ps1LoopSmooth) ;

				static StoreM TableM ;

				if (ImportM == false) ResumM (PathToOutput, paramsP11, ExportM, &TableM) ;

				//ExportPowerSpectraNoResum (PathToOutput, 0, &PsLinear) ;
				//ExportPowerSpectraNoResum (PathToOutput, 1, &Ps1Loop) ;
				
				ResumSmoothOscPowerSpectra (PathToOutput, paramsP11, &PsLinear, &Ps1Loop, &PsLinearSmooth, &Ps1LoopSmooth, &TableM) ;
				//ResumPowerSpectra (PathToOutput, paramsP11, &PsLinear, &Ps1Loop, ImportM, ExportM, &TableM) ;
			}

			else {
				ExportPowerSpectraNoResum (PathToOutput, 0, &PsLinear) ;
				ExportPowerSpectraNoResum (PathToOutput, 1, &Ps1Loop) ;

			}
		}

		if (ComputeBispectrum == true) {
			ComputeBispectrumMonopole (PathToOutput, PathToTriangles, paramsP11, nbar, aperp, apar) ;
		}
	}

	return 0 ;

}