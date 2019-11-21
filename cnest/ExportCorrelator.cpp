#include "CBiRd.h"

void ExportCorrelator (const string & Path, const char * Name, const size_t & Nmax, Coordinates * x, Correlator * O, const YesNo & ShotNoise, const unsigned int & Nlmax) {

	int pre = 12 ; // white spacing

	ostringstream filename ;
	filename << Path << "/" << Name << ".dat" ;

	ofstream write (filename.str(), ios::out | ios::trunc) ;

	write << "#" << setw(pre-1) << "k" << setw(pre) 
		<< "0-1" << setw(pre) << "0-b1" << setw(pre) << "0-b1*b1" << setw(pre) 
		<< "1" << setw(pre) << "b1"  << setw(pre) << "b2"  << setw(pre) << "b3"  << setw(pre) << "b4"  << setw(pre)
		<< "b1*b1"  << setw(pre) << "b1*b2"  << setw(pre) << "b1*b3"  << setw(pre) << "b1*b4"  << setw(pre)
		<< "b2*b2"  << setw(pre) << "b2*b4"  << setw(pre) << "b4*b4" << setw(pre) 
		<< "b1*b5"  << setw(pre) << "b1*b6"  << setw(pre) << "b1*b7" << setw(pre) 
		<< "b5"  << setw(pre) << "b6"  << setw(pre) << "b7"  << setw(pre) << endl ;
		//<< "b8"  << setw(pre) << "b9"  << setw(pre) << "b10" << endl ;

	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		for (unsigned int m = 0 ; m < Nmax ; m++) {
			write << setw(pre) << (*x)[m] << " " ;
			for (unsigned int n = 0 ; n < Np ; n++) {
				if (ShotNoise == true && n >= 3 && n < 15) write << setw(pre-1) << (*O)[i][n][m] - (*O)[i][n][0] << " " ;
				else write << setw(pre-1) << (*O)[i][n][m] << " " ;
			}
			write << endl ;
		}
	}
	write.close() ;
}
