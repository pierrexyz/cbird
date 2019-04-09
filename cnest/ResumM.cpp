#include "fftlog.h"
#include "ResumEFT.h"
#include "kout.h"


void ResumM (const string & PathToFolder, const ParamsP11 & params, const bool & ExportM, StoreM * TableM) {

	// FFTLOG
	double kpmin = 0.001, kpmax = 1000. ; 	// Changing kpmin from 0.001 to 0.1 makes a 0.2% difference
											// I don't know which one is the best, but I keep 0.001 for maximal periodic boundaries
	double dl = (log(kpmax)-log(kpmin))/(Nkp-1.) ;
	double kplist[Nkp] ;
	for (unsigned int i = 0 ; i< Nkp ; i++) kplist[i] = kpmin*exp(dl*i) ;

	double L = log(kplist[Nkp-1]/kplist[0]) * Nkp/(Nkp-1.) ;
	

	double X1Infinity ;
	InterpFunc InterpX1, InterpY1 ;
	//GetX1Y1(params, InterpX1, InterpY1, X1Infinity) ;

	double X1Table[Nkp], Y1Table[Nkp] ;

	for (unsigned int m = 0 ; m < Nkp ; m++) {
		X1Table[m] = ResumX1(kplist[m],params) ;
		Y1Table[m] = ResumY1(kplist[m],params) ;
	}

	InterpX1.accel = gsl_interp_accel_alloc () ;
	InterpX1.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
	gsl_spline_init (InterpX1.interp, kplist, X1Table, Nkp) ;

	InterpY1.accel = gsl_interp_accel_alloc () ;
	InterpY1.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
	gsl_spline_init (InterpY1.interp, kplist, Y1Table, Nkp) ;

	X1Infinity = ResumX1(qInfinity,params) ;

	// Compute of M^o_l,lp (k,kp) ; o = 0,1 ; l = 0,2,4 ; lp = 0,2,4,6,8 ; k given by kout.h
	for (unsigned int order = 0 ; order < 2 ; order++) {

		size_t Nlp ;
		if (order == 0) Nlp = 3 ; // Nlp = 5 : sub-percent difference
		if (order == 1) Nlp = 3 ;
		
		for (unsigned int j = 0 ; j < Nlp ; j++) { // loop over lp

			int lp = 2*j ;

			//low-ringing condition
			double kcrc = goodkr(Nkp, lp+0.5, 0., L, 1.) ;
			dcomplex u[Nkp] ;
			compute_u_coefficients(Nkp, lp+0.5, 0., L, kcrc, u) ;

			for (unsigned int i = 0 ; i < Nlout ; i++) { // loop over l

				int l = 2*i ;

				for (unsigned int m = 0 ; m < Nout ; m++) {	// loop over k

					double k = kout[m] ;

					dcomplex Mhat[Nkp] ;

					for (unsigned int n = 0 ; n < Nkp ; n++) {

						double ChangeToSphericalBessel = pow(kplist[n],-0.5) ;
						

						if (order == 0) Mhat[n] = ChangeToSphericalBessel *  ResumQ0(l,lp,k,kplist[n],InterpX1,InterpY1,X1Infinity,params.f) ;
						if (order == 1) Mhat[n] = ChangeToSphericalBessel *  ResumQ1(l,lp,k,kplist[n],InterpX1,InterpY1,X1Infinity,params.f) ;
					}


    				double kres[Nkp] ;
    				dcomplex result[Nkp] ;
    				
    				fht(Nkp, kplist, Mhat, kres, result, lp+0.5, 0., kcrc, false, u) ;	
    				//fht(Nkp, kplist, Mhat, kres, result, lp+0.5, 0.) ;	

    				for (unsigned int n = 0 ; n < Nkp ; n++) {
    					(*TableM)[0][order][i][j][m][n] = kres[n] ;
    					(*TableM)[1][order][i][j][m][n] = sqrt(M_PI/2.)* pow(kres[n], -1.5)* real(result[n]) ;
    				}

    				if (ExportM == true) {
    					ostringstream filename ;
    					filename << PathToFolder << setprecision(3) << "resum_data/M" << order << "_" << k << "_" << l << "_" << lp << ".dat" ;
    					ofstream write (filename.str(), ios::out | ios::trunc) ;
    					for (unsigned int n = 0 ; n < Nkp ; n++) write << setw(12) << kres[n] << " " << setw(12) << (*TableM)[1][order][i][j][m][n] << endl ;
    					write.close() ;

    					ostringstream filename2 ;
    					filename2 << PathToFolder << setprecision(3) << "resum_data/QInfinity" << order << "_" << k << "_" << l << "_" << lp << ".dat" ;
    					ofstream write2 (filename2.str(), ios::out | ios::trunc) ;
    					write2 << QInfinity(order,k,l,lp,X1Infinity,params.f) << endl ;
    					write2.close() ;
    				}

				}
			}
		}

	}

	UnloadInterp(InterpX1) ;
	UnloadInterp(InterpY1) ;

}



double QInfinity (const unsigned & Morder, const double & k, const unsigned & l, const unsigned & lp, const double & X1qInfinity, const double & f1) {
	if (Morder == 1) return 0.5*(2*l+1.) * exp(-0.5*k*k*X1qInfinity) * ( (1.+0.5*k*k*X1qInfinity)*(*ResumI0[l/2][lp/2])(k,X1qInfinity,f1) + 0.5*k*k*f1*(2.+f1)*X1qInfinity*(*ResumI2[l/2][lp/2])(k,X1qInfinity,f1) )  ;
	if (Morder == 0) return 0.5*(2*l+1.) * exp(-0.5*k*k*X1qInfinity) * (*ResumI0[l/2][lp/2])(k,X1qInfinity,f1) ;
}

double LoadQInfinity (const string & PathToFolder, const unsigned & Morder, const double & k, const unsigned & l, const unsigned & lp) {

	ostringstream filename ;
    filename << setprecision(3) << PathToFolder << "resum_data/QInfinity" << Morder << "_" << k << "_" << l << "_" << lp << ".dat" ;
	ifstream data(filename.str(), ios::in) ;

	if (!data.is_open()) { 
		cerr << "Problem loading file:" << filename.str() << endl ; 
		exit(EXIT_FAILURE) ; 
	}

	else {
		double Qinf ;
		data >> Qinf ;
		data.close() ;
		return Qinf ;
	}
}