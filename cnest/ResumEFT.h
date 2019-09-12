#ifndef RESUMEFT_H
#define RESUMEFT_H

#include "CBiRd.h"
#include "kout.h"


struct ParamsResumIntegr {
	InterpFunc InterpP, InterpM ;
} ;



//////////////////////////////////////////////////////////////////
// Functions needed to create the M-matrices.
//////////////////////////////////////////////////////////////////

// ResumI.cpp
typedef double ResumI (const double & k, const double & X1q, const double & f1) ;

ResumI I00, I02, I04, I06, I08, I22, I24, I26, I28, I44, I46, I48, I66, I68, I88 ;
ResumI J00, J02, J04, J22, J24, J44 ;

ResumI * const ResumI0[5][5] = { 
	{ I00, I02, I04, I06, I08 } ,
	{ I02, I22, I24, I26, I28 } ,
	{ I04, I24, I44, I46, I48 } ,
	{ I06, I26, I46, I66, I68 } ,
	{ I08, I28, I48, I68, I88 } 
} ;

ResumI * const ResumI2[3][3] = {
	{ J00, J02, J04 } ,
	{ J02, J22, J24 } ,
	{ J04, J24, J44 } ,
} ; 

typedef complex<double> dcomplex ;

const double qInfinity = 10000. ;
const double alpha = 1./3. ;
const double beta = 1./3. ;


const size_t Nlout = 3 ; // l = 0,2,4

// M^o_l,lp (k,kp) : StoreM [new kp:0 | M:1][o][l][lp][k][kp]
typedef double StoreM [2][2][Nlout][Nl][Nout][Nkp] ;

// ResumX1Y1.cpp
double ResumX1 (const double & q, const ParamsP11 & InterpP11) ;
double ResumY1 (const double & q, const ParamsP11 & InterpP11) ;
void GetX1Y1 (const ParamsP11 & InterpP11, InterpFunc & InterpX1, InterpFunc & InterpY1, double & X1Infinity) ;

// ResumQ.cpp
double ResumQ0 (const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & X1Infinity, const double & f1) ;
double ResumQ1(const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & X1Infinity, const double & f1) ;

// ResumM.cpp
void ResumM (const string & PathToFolder, const ParamsP11 & params, const bool & ExportM, StoreM * TableM) ;
double QInfinity (const unsigned & order, const double & k, const unsigned & l, const unsigned & lp, const double & X1Infinity, const double & f1) ;
double LoadQInfinity (const string & PathToFolder, const unsigned & Morder, const double & k, const unsigned & l, const unsigned & lp) ;


//////////////////////////////////////////////////////////////////
// IR-resummation
//////////////////////////////////////////////////////////////////

// LoadResummation.cpp
void LoadM (const string & PathToFolder, const unsigned int & order, const double & k, const unsigned int & l, const unsigned int & lp, InterpFunc & InterpM) ;
void UnloadInterp (InterpFunc & params) ;

// Resummation.cpp
double P (const double & q, const InterpFunc & InterpP) ;
double M (const double & q, const InterpFunc & InterpM) ;
double X1 (const double & q, const InterpFunc & InterpX1) ;
double Y1 (const double & q, const InterpFunc & InterpX1) ;

void ResumPowerSpectra (const string & PathToFolder, const ParamsP11 & InterpP11, PowerSpectraNoResum * Linear, PowerSpectraNoResum * Loop, const YesNo & ImportM, const YesNo & ExportM, StoreM * TableM) ;
double Resum_Integrand_GSL (double q, void * params) ;

// ResummationSmoothOsc.cpp
void ResumSmoothOscPowerSpectra (const string & PathToFolder, const ParamsP11 & params, PowerSpectraNoResum * Linear, PowerSpectraNoResum * Loop, PowerSpectraNoResum * LinearSmooth, PowerSpectraNoResum * LoopSmooth, StoreM * TableM) ;

#endif
