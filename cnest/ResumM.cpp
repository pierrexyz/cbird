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
						// ostringstream filename ;
						// filename << PathToFolder << setprecision(3) << "/resum_data/M" << order << "_" << k << "_" << l << "_" << lp << ".dat" ;
						// ofstream write (filename.str(), ios::out | ios::trunc) ;
						// for (unsigned int n = 0 ; n < Nkp ; n++) write << setw(12) << kres[n] << " " << setw(12) << (*TableM)[1][order][i][j][m][n] << endl ;
						// write.close() ;

						// ostringstream filename2 ;
    					// filename2 << PathToFolder << setprecision(3) << "/resum_data/Q0_" << k << "_" << l << "_" << lp << ".dat" ;
    					// ofstream write2 (filename2.str(), ios::out | ios::trunc) ;
    					// for (unsigned int n = 0 ; n < Nkp ; n++) write2 << setw(12) << kplist[n] << " " << setw(12) << ResumQ0(l,lp,k,kplist[n],InterpX1,InterpY1,X1Infinity,params.f) << endl ;
    					// write2.close() ;

						// ostringstream filename3 ;
    					// filename3 << PathToFolder << setprecision(3) << "/resum_data/Q1_" << k << "_" << l << "_" << lp << ".dat" ;
    					// ofstream write3 (filename3.str(), ios::out | ios::trunc) ;
    					// for (unsigned int n = 0 ; n < Nkp ; n++) write << setw(12) << kplist[n] << " " << setw(12) << ResumQ1(l,lp,k,kplist[n],InterpX1,InterpY1,X1Infinity,params.f) << endl ;
    					// write3.close() ;



    					// ostringstream filename2 ;
    					// filename2 << PathToFolder << setprecision(3) << "resum_data/QInfinity" << order << "_" << k << "_" << l << "_" << lp << ".dat" ;
    					// ofstream write2 (filename2.str(), ios::out | ios::trunc) ;
    					// write2 << QInfinity(order,k,l,lp,X1Infinity,params.f) << endl ;
    					// write2.close() ;
    				}

				}
			}
		}

	}

	UnloadInterp(InterpX1) ;
	UnloadInterp(InterpY1) ;

}



double QInfinity (const unsigned & Morder, const double & k, const unsigned & l, const unsigned & lp, const double & X1Infinity, const double & f1) {

	int i = l/2, j = lp/2 ;

	if (Morder == 1) {
		
switch (i) { 
	case 0:
		switch (j) { 
			case 0:
				return (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*(2*f1*(2 + f1)*X1Infinity*pow(k,2)*((-6*f1*(2 + f1)*X1Infinity*pow(k,2)*(8 + beta*f1*(2 + f1)*X1Infinity*pow(k,2)*(4 + beta*f1*(2 + f1)*X1Infinity*pow(k,2))))/5. + (2*(48 + beta*f1*(2 + f1)*X1Infinity*pow(k,2)*(24 + beta*f1*(2 + f1)*X1Infinity*pow(k,2)*(6 + beta*f1*(2 + f1)*X1Infinity*pow(k,2)))))/3. + (6*pow(f1,2)*pow(2 + f1,2)*(2 + beta*f1*(2 + f1)*X1Infinity*pow(k,2))*pow(k,4)*pow(X1Infinity,2))/7. - (2*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1Infinity,3))/9.) + 4*(1 + (X1Infinity*pow(k,2))/2.)*(96 + 16*(-1 + 3*beta)*f1*(2 + f1)*X1Infinity*pow(k,2) + (4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 + f1,2)*pow(k,4)*pow(X1Infinity,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + beta)*beta))*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1Infinity,3))/35.)))/384.  ;
				break ;

			case 1:
				return  (f1*(2 + f1)*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,4)*pow(X1Infinity,2)*(-5544 - 792*f1*(6 + (-3 + 7*beta)*X1Infinity*pow(k,2)) + 4*X1Infinity*pow(f1,3)*pow(k,2)*(660 - 445*X1Infinity*pow(k,2) - 1881*X1Infinity*pow(beta,2)*pow(k,2) + 462*X1Infinity*pow(beta,3)*pow(k,2) + 396*beta*(-3 + 4*X1Infinity*pow(k,2))) + 3*X1Infinity*pow(f1,4)*pow(k,2)*(220 - 615*X1Infinity*pow(k,2) - 2607*X1Infinity*pow(beta,2)*pow(k,2) + 924*X1Infinity*pow(beta,3)*pow(k,2) + 198*beta*(-2 + 11*X1Infinity*pow(k,2))) + 6*(-140 + 495*beta - 594*pow(beta,2) + 231*pow(beta,3))*pow(f1,5)*pow(k,4)*pow(X1Infinity,2) + (-140 + 495*beta - 594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,4)*pow(X1Infinity,2) - 132*pow(f1,2)*(18 + (-29 + 57*beta)*X1Infinity*pow(k,2) + (5 - 18*beta + 21*pow(beta,2))*pow(k,4)*pow(X1Infinity,2))))/166320.  ;
				break ;

			case 2:
				return  -(pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,2)*pow(2 + f1,2)*pow(k,4)*pow(X1Infinity,2)*(572 + 26*X1Infinity*(-11 + (-20 + 22*beta)*f1 + (-10 + 11*beta)*pow(f1,2))*pow(k,2) + f1*(2 + f1)*(65 + 140*f1 + 143*f1*(2 + f1)*pow(beta,2) + 70*pow(f1,2) - 13*beta*(11 + 30*f1 + 15*pow(f1,2)))*pow(k,4)*pow(X1Infinity,2)))/180180. ;
				break ;

		}
		break ;

	case 1:
		switch (j) { 
			case 0:
				return (f1*(2 + f1)*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,4)*pow(X1Infinity,2)*(-5544 - 792*f1*(6 + (-3 + 7*beta)*X1Infinity*pow(k,2)) + 4*X1Infinity*pow(f1,3)*pow(k,2)*(660 - 445*X1Infinity*pow(k,2) - 1881*X1Infinity*pow(beta,2)*pow(k,2) + 462*X1Infinity*pow(beta,3)*pow(k,2) + 396*beta*(-3 + 4*X1Infinity*pow(k,2))) + 3*X1Infinity*pow(f1,4)*pow(k,2)*(220 - 615*X1Infinity*pow(k,2) - 2607*X1Infinity*pow(beta,2)*pow(k,2) + 924*X1Infinity*pow(beta,3)*pow(k,2) + 198*beta*(-2 + 11*X1Infinity*pow(k,2))) + 6*(-140 + 495*beta - 594*pow(beta,2) + 231*pow(beta,3))*pow(f1,5)*pow(k,4)*pow(X1Infinity,2) + (-140 + 495*beta - 594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,4)*pow(X1Infinity,2) - 132*pow(f1,2)*(18 + (-29 + 57*beta)*X1Infinity*pow(k,2) + (5 - 18*beta + 21*pow(beta,2))*pow(k,4)*pow(X1Infinity,2))))/33264.  ;
				break ;

			case 1:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*(864864 + 432432*(1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2) + 5148*f1*(2 + f1)*(-22 + 42*beta - 18*f1 + 21*f1*(2 + f1)*pow(beta,2) - 9*pow(f1,2))*pow(k,4)*pow(X1Infinity,2) + 78*pow(f1,2)*(297 + 340*f1 + 693*pow(beta,2) + 231*f1*(2 + f1)*pow(beta,3) + 170*pow(f1,2) - 33*beta*(22 + 18*f1 + 9*pow(f1,2)))*pow(2 + f1,2)*pow(k,6)*pow(X1Infinity,3) + (-1287*pow(beta,2)*(11 + 18*f1 + 9*pow(f1,2)) + 429*pow(beta,3)*(21 + 22*f1 + 11*pow(f1,2)) + 117*beta*(99 + 170*f1 + 85*pow(f1,2)) - 5*(663 + 1162*f1 + 581*pow(f1,2)))*pow(f1,3)*pow(2 + f1,3)*pow(k,8)*pow(X1Infinity,4)))/864864.  ;
				break ;

			case 2:
				return  (f1*(2 + f1)*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,4)*pow(X1Infinity,2)*(-3432 - 104*f1*(34 + (-17 + 33*beta)*X1Infinity*pow(k,2)) + 4*X1Infinity*pow(f1,3)*pow(k,2)*(580 - 425*X1Infinity*pow(k,2) - 1313*X1Infinity*pow(beta,2)*pow(k,2) + 286*X1Infinity*pow(beta,3)*pow(k,2) + 4*beta*(-221 + 328*X1Infinity*pow(k,2))) + X1Infinity*pow(f1,4)*pow(k,2)*(580 - 1825*X1Infinity*pow(k,2) - 5733*X1Infinity*pow(beta,2)*pow(k,2) + 1716*X1Infinity*pow(beta,3)*pow(k,2) + beta*(-884 + 5662*X1Infinity*pow(k,2))) + 6*(-140 + 435*beta - 442*pow(beta,2) + 143*pow(beta,3))*pow(f1,5)*pow(k,4)*pow(X1Infinity,2) + (-140 + 435*beta - 442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,4)*pow(X1Infinity,2) - 4*pow(f1,2)*(442 + (-801 + 1313*beta)*X1Infinity*pow(k,2) + (145 - 442*beta + 429*pow(beta,2))*pow(k,4)*pow(X1Infinity,2))))/72072. ;
				break ;

		}
		break ;

	case 2:
		switch (j) { 
			case 0:
				return -(pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,2)*pow(2 + f1,2)*pow(k,4)*pow(X1Infinity,2)*(572 + 26*X1Infinity*(-11 + (-20 + 22*beta)*f1 + (-10 + 11*beta)*pow(f1,2))*pow(k,2) + f1*(2 + f1)*(65 + 140*f1 + 143*f1*(2 + f1)*pow(beta,2) + 70*pow(f1,2) - 13*beta*(11 + 30*f1 + 15*pow(f1,2)))*pow(k,4)*pow(X1Infinity,2)))/20020.  ;
				break ;

			case 1:
				return  (f1*(2 + f1)*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,4)*pow(X1Infinity,2)*(-3432 - 104*f1*(34 + (-17 + 33*beta)*X1Infinity*pow(k,2)) + 4*X1Infinity*pow(f1,3)*pow(k,2)*(580 - 425*X1Infinity*pow(k,2) - 1313*X1Infinity*pow(beta,2)*pow(k,2) + 286*X1Infinity*pow(beta,3)*pow(k,2) + 4*beta*(-221 + 328*X1Infinity*pow(k,2))) + X1Infinity*pow(f1,4)*pow(k,2)*(580 - 1825*X1Infinity*pow(k,2) - 5733*X1Infinity*pow(beta,2)*pow(k,2) + 1716*X1Infinity*pow(beta,3)*pow(k,2) + beta*(-884 + 5662*X1Infinity*pow(k,2))) + 6*(-140 + 435*beta - 442*pow(beta,2) + 143*pow(beta,3))*pow(f1,5)*pow(k,4)*pow(X1Infinity,2) + (-140 + 435*beta - 442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,4)*pow(X1Infinity,2) - 4*pow(f1,2)*(442 + (-801 + 1313*beta)*X1Infinity*pow(k,2) + (145 - 442*beta + 429*pow(beta,2))*pow(k,4)*pow(X1Infinity,2))))/40040.  ;
				break ;

			case 2:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*(8168160 + 4084080*(1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2) + 204*f1*(2 + f1)*(10010*beta + 5005*f1*(2 + f1)*pow(beta,2) - 3*(1690 + 1286*f1 + 643*pow(f1,2)))*pow(k,4)*pow(X1Infinity,2) + 34*pow(f1,2)*(5787 + 6540*f1 + 15015*pow(beta,2) + 5005*f1*(2 + f1)*pow(beta,3) + 3270*pow(f1,2) - 9*beta*(1690 + 1286*f1 + 643*pow(f1,2)))*pow(2 + f1,2)*pow(k,6)*pow(X1Infinity,3) + (1105*pow(beta,3)*(77 + 78*f1 + 39*pow(f1,2)) + 153*beta*(643 + 1090*f1 + 545*pow(f1,2)) - 153*pow(beta,2)*(845 + 1286*f1 + 643*pow(f1,2)) - 15*(1853 + 3318*f1 + 1659*pow(f1,2)))*pow(f1,3)*pow(2 + f1,3)*pow(k,8)*pow(X1Infinity,4)))/8.16816e6 ;
				break ;

		}
		break ;

}



//return 0.5*(2*l+1.) * exp(-0.5*k*k*X1Infinity) * ( (1.+0.5*k*k*X1Infinity)*(*ResumI0[l/2][lp/2])(k,X1Infinity,f1) + 0.5*k*k*f1*(2.+f1)*X1Infinity*(*ResumI2[l/2][lp/2])(k,X1Infinity,f1) )  ;
	}

	if (Morder == 0) {

switch (i) { 
	case 0:
		switch (j) { 
			case 0:
				return (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*(96 + 16*(-1 + 3*beta)*f1*(2 + f1)*X1Infinity*pow(k,2) + (4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 + f1,2)*pow(k,4)*pow(X1Infinity,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + beta)*beta))*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1Infinity,3))/35.))/96.  ;
				break ;

			case 1:
				return  (f1*(2 + f1)*X1Infinity*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,2)*(-168 + f1*(2 + f1)*X1Infinity*pow(k,2)*(36 - 84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + f1)*X1Infinity*pow(k,2))))/2520.  ;
				break ;

			case 2:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,2)*pow(2 + f1,2)*(22 + (-5 + 11*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))*pow(k,4)*pow(X1Infinity,2))/6930.  ;
				break ;

			case 3:
				return  -(pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1Infinity,3))/9009.  ;
				break ;

			case 4:
				return  0 ;
				break ;

		}
		break ;

	case 1:
		switch (j) { 
			case 0:
				return (f1*(2 + f1)*X1Infinity*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,2)*(-168 + f1*(2 + f1)*X1Infinity*pow(k,2)*(36 - 84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + f1)*X1Infinity*pow(k,2))))/504.  ;
				break ;

			case 1:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*(11088 + 528*(-11 + 21*beta)*f1*X1Infinity*pow(k,2) + 264*X1Infinity*pow(f1,2)*pow(k,2)*(-11 + 21*beta + X1Infinity*(9 - 22*beta + 21*pow(beta,2))*pow(k,2)) + 8*pow(f1,3)*(33*(9 - 22*beta + 21*pow(beta,2)) + X1Infinity*(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2))*pow(k,4)*pow(X1Infinity,2) + 6*pow(f1,4)*(11*(9 - 22*beta + 21*pow(beta,2)) + 2*X1Infinity*(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2))*pow(k,4)*pow(X1Infinity,2) + 6*(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3)))/11088.  ;
				break ;

			case 2:
				return  (f1*(2 + f1)*X1Infinity*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,2)*(-3432 + f1*(2 + f1)*X1Infinity*pow(k,2)*(884 - 1716*beta + (-145 + 13*(34 - 33*beta)*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))))/36036.  ;
				break ;

			case 3:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,2)*pow(2 + f1,2)*(90 + (-23 + 45*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))*pow(k,4)*pow(X1Infinity,2))/18018.  ;
				break ;

			case 4:
				return  (-4*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1Infinity,3))/21879. ;
				break ;

		}
		break ;

	case 2:
		switch (j) { 
			case 0:
				return (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,2)*pow(2 + f1,2)*(22 + (-5 + 11*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))*pow(k,4)*pow(X1Infinity,2))/770.  ;
				break ;

			case 1:
				return  (f1*(2 + f1)*X1Infinity*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,2)*(-3432 + f1*(2 + f1)*X1Infinity*pow(k,2)*(884 - 1716*beta + (-145 + 13*(34 - 33*beta)*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))))/20020.  ;
				break ;

			case 2:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*(240240 + 3120*(-39 + 77*beta)*f1*X1Infinity*pow(k,2) + 24*X1Infinity*pow(f1,2)*pow(k,2)*(65*(-39 + 77*beta) + X1Infinity*(1929 - 5070*beta + 5005*pow(beta,2))*pow(k,2)) + 8*pow(f1,3)*(3*(1929 - 5070*beta + 5005*pow(beta,2)) + X1Infinity*(-1635 + 5787*beta - 7605*pow(beta,2) + 5005*pow(beta,3))*pow(k,2))*pow(k,4)*pow(X1Infinity,2) + 6*pow(f1,4)*(1929 - 5070*beta + 5005*pow(beta,2) + 2*X1Infinity*(-1635 + 5787*beta - 7605*pow(beta,2) + 5005*pow(beta,3))*pow(k,2))*pow(k,4)*pow(X1Infinity,2) + 6*(-1635 + 5787*beta - 7605*pow(beta,2) + 5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-1635 + 5787*beta - 7605*pow(beta,2) + 5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3)))/240240.  ;
				break ;

			case 3:
				return  (f1*(2 + f1)*X1Infinity*pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(k,2)*(-14280 + f1*(2 + f1)*X1Infinity*pow(k,2)*(3604 - 7140*beta + (-569 + 17*(106 - 105*beta)*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))))/136136.  ;
				break ;

			case 4:
				return  (pow(E,-((1 + beta*f1*(2 + f1))*X1Infinity*pow(k,2))/2.)*pow(f1,2)*pow(2 + f1,2)*(266 + (-67 + 133*beta)*f1*(2 + f1)*X1Infinity*pow(k,2))*pow(k,4)*pow(X1Infinity,2))/46189. ;
				break ;

		}
		break ;

}


		//return 0.5*(2*l+1.) * exp(-0.5*k*k*X1Infinity) * (*ResumI0[l/2][lp/2])(k,X1Infinity,f1) ;
	}
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