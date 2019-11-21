#include "CBiRd.h"

static YesNo GetAnswer (const string & Answer) {
	if (Answer == "Y" || Answer == "y" || Answer == "yes" || Answer == "Yes" || Answer == "YES") return 1 ;
	if (Answer == "N" || Answer == "n" || Answer == "no" || Answer == "No" || Answer == "NO") return 0 ;
	else return 1 ; 
}

void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, 
	string & PathToOutput, string & PathToLinearPowerSpectrum, 
	YesNo & ExportPowerSpectrumNoResum, YesNo & ExportCorrFuncNoResum, YesNo & ResumPowerSpectrum,
	double & aperp, double & apar, unsigned int & Nlout) {
///// Reading ini file /////
	ifstream inFile (ConfigFile) ;    
	char line[1024] ;
	char delim[] = " =" ;
	char * token ;
	
	if (!inFile.is_open()) {
		cout << "Could not open the input file: " << ConfigFile << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	while (inFile) {
		
		inFile.getline(line, 1024) ;
		
		if (strlen(line) == 0 or line[0]=='#' or line[0] == '[') continue ;
		
		token = strtok(line, delim) ;
		
		if (strcmp(token, "PathToOutput") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToOutput = token ;
		else if (strcmp(token, "PathToLinearPowerSpectrum") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToLinearPowerSpectrum = token ;

		else if (strcmp(token, "knl") == 0) while ((token = strtok(NULL, delim)) != NULL) knl = atof(token) ;
		else if (strcmp(token, "km") == 0) while ((token = strtok(NULL, delim)) != NULL) km = atof(token) ;
		else if (strcmp(token, "nbar") == 0) while ((token = strtok(NULL, delim)) != NULL) nbar = atof(token) ;

		else if (strcmp(token, "ExportPowerSpectrumNoResum") == 0) while ((token = strtok(NULL, delim)) != NULL) ExportPowerSpectrumNoResum = GetAnswer(token) ;
		else if (strcmp(token, "ExportCorrFuncNoResum") == 0) while ((token = strtok(NULL, delim)) != NULL) ExportCorrFuncNoResum = GetAnswer(token) ;
		else if (strcmp(token, "ResumPowerSpectrum") == 0) while ((token = strtok(NULL, delim)) != NULL) ResumPowerSpectrum = GetAnswer(token) ;

		else if (strcmp(token, "z_pk") == 0) while ((token = strtok(NULL, delim)) != NULL) z0 = atof(token) ;
		else if (strcmp(token, "ln10^{10}A_s") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[0] = atof(token) ;
		else if (strcmp(token, "A_s") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[0] = log( atof(token) *1e10 ) ;
		else if (strcmp(token, "n_s") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[1] = atof(token) ;
		else if (strcmp(token, "h") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[2] = atof(token) ;
		else if (strcmp(token, "omega_b") == 0) while ((token = strtok(NULL, delim)) != NULL) { cosmo[3] = atof(token) ; }
		else if (strcmp(token, "omega_cdm") == 0) while ((token = strtok(NULL, delim)) != NULL) { cosmo[4] = atof(token) ; }
		else if (strcmp(token, "N_ncdm") == 0) while ((token = strtok(NULL, delim)) != NULL) { cosmo[5] = atof(token) ; }
		else if (strcmp(token, "Sum_mnu") == 0) while ((token = strtok(NULL, delim)) != NULL) { cosmo[6] = atof(token) ; }
		
		else if (strcmp(token, "aperp") == 0) while ((token = strtok(NULL, delim)) != NULL) aperp = atof(token) ;
		else if (strcmp(token, "apar") == 0) while ((token = strtok(NULL, delim)) != NULL) apar = atof(token) ;

		// else if (strcmp(token, "Multipoles") == 0) while ((token = strtok(NULL, delim)) != NULL) Nlout = atof(token) ;
		}
	
	inFile.close() ;
}

