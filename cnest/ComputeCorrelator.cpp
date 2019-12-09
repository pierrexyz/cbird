#include "CBiRd.h"

static dcomplex Mat22[Nl][N22][(NFFT+1)*(NFFT+1)] ;
static dcomplex Mat13[Nl][N13][(NFFT+1)*(NFFT+1)] ;

// 'M'atrix for 'P'ower-spectrum to 'C'orrelation-function Hankel-transform
// Ps is a power-law 'n1' cosmology of multipole 'l'
dcomplex MPC (const unsigned int & l, const dcomplex & n1) {  
	return pow(2., -2.*n1) * pow(M_PI, -1.5) * MyGamma(1.5 + l/2. - n1) / MyGamma(l/2. + n1) ;
}

void ComputeCorrelator (const ParamsP11 & p, const ParamsP11 & psmooth, const size_t & Nq, Coordinates * s, Correlator * Cf , Correlator * Cfsmooth , const size_t & Nk, Coordinates * k , Correlator * Ps, const unsigned int & Nlmax) {
	
	// FFTLog decomposition
	dcomplex Coef[NFFT+1], Pow[NFFT+1] ;
	CoefPow(Coef, Pow, p) ;

	Eigen::Map<Eigen::ArrayXcd> ArrPow(Pow, NFFT+1) ;
	Eigen::Map<Eigen::ArrayXcd> ArrCoef(Coef, NFFT+1) ;

	/****** Power spectrum ******/

	// Terms c_m k^p_m in the sum of P11 (from FFTLog decomposition)
	Eigen::ArrayXcd ArrTermsPs[Nk] ;
	for (unsigned int m = 0 ; m < Nk ; m++) {
		ArrTermsPs[m] = ArrPow ;
		ArrTermsPs[m] *= log((*k)[m]) ;
		ArrTermsPs[m] = ArrTermsPs[m].exp() ;
		ArrTermsPs[m] *= ArrCoef ;
	}

	/* Ps-linear */
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		for (unsigned int m = 0 ; m < Nk; m++) {
			double P11k = P11((*k)[m], p) ;
			(*Ps)[i][0][m] = m4[i] *p.f*p.f *P11k ; // 1
			(*Ps)[i][1][m] = m2[i] *2.*p.f *P11k ; 	// b1
			(*Ps)[i][2][m] = m0[i] *P11k ; 			// b1*b1
		}
	}

	/* Ps-counterterm */
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		for (unsigned int m = 0 ; m < Nk ; m++) {
			double P11k = P11((*k)[m], p) ;
			double k2 = pow((*k)[m], 2) ;
			(*Ps)[i][15][m] = m0[i] * 2. *k2 *P11k ;			// b1*b5
			(*Ps)[i][16][m] = m2[i] * 2. *k2 *P11k ;			// b1*b6
			(*Ps)[i][17][m] = m4[i] * 2. *k2 *P11k ; 			// b1*b7
			(*Ps)[i][18][m] = m2[i] * 2. *p.f *k2 *P11k ;		// b5
			(*Ps)[i][19][m] = m4[i] * 2. *p.f *k2 *P11k ;		// b6
			(*Ps)[i][20][m] = m6[i] * 2. *p.f *k2 *P11k ;		// b7
		}
	}

	/* Ps-1loop */
	double P22[N22][Nk], P13[N13][Nk] ;

	// M13 from Loop
	for (unsigned int j = 0 ; j < N13 ; j++)
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
				Mat13[0][j][m1] = M13(j, -0.5*Pow[m1], p.f) ; 
	// P13
	for (unsigned int j = 0 ; j < N13 ; j++) {
		Eigen::RowVectorXcd M13p(NFFT+1) ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) M13p(m1) = Mat13[0][j][m1] ;
		for (unsigned int i = 0 ; i < Nk ; i++) {
			dcomplex tmp = M13p*ArrTermsPs[i].matrix() ;
			P13[j][i] = pow((*k)[i],3)*P11((*k)[i],p)*tmp.real() ;
		}
	}
	
	// Auxiliary function involving gamma functions for M22
	dcomplex In1n2[NFFT+1][NFFT+1] ;
	for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
		for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
			In1n2[m1][m2] = I(-0.5*Pow[m1], -0.5*Pow[m2]) ;
	// M22 from Loop
	for (unsigned int j = 0 ; j < N22 ; j++)
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
			for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
				Mat22[0][j][m1+(NFFT+1)*m2] = M22(j, -0.5*Pow[m1], -0.5*Pow[m2], In1n2[m1][m2], p.f) ; 
	// Shot-noise in P22 to subtract
	// Eigen::ArrayXcd ArrTerm0 ;
	// ArrTerm0 = ArrPow ;
	// ArrTerm0 *= log(k0) ;
	// ArrTerm0 = ArrTerm0.exp() ;
	// ArrTerm0 *= ArrCoef ;
	// P22
	for (unsigned int j = 0 ; j < N22 ; j++) {
		Eigen::MatrixXcd M22p(NFFT+1,NFFT+1) ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
			for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
				M22p(m1,m2) = Mat22[0][j][m1+(NFFT+1)*m2] ;

		// dcomplex tmp0 = ArrTerm0.matrix().transpose() * M22p * ArrTerm0.matrix() ;
		// double ShotNoise = pow(k0,3)*tmp0.real() ;	

		for (unsigned int i = 0 ; i < Nk ; i++) {
			dcomplex tmp = ArrTermsPs[i].matrix().transpose() * M22p * ArrTermsPs[i].matrix() ;
			P22[j][i] = pow((*k)[i],3)*tmp.real() ;
			// P22[j][i] -= ShotNoise ;
		}
	}

	// Ps_l
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		for (unsigned int m = 0 ; m < Nk ; m++) {
			(*Ps)[i][3][m] = (*Multipole22[17])[i]*P22[17][m] + (*Multipole22[18])[i]*P22[18][m] + (*Multipole22[19])[i]*P22[19][m] + (*Multipole13[6])[i]*P13[6][m] + (*Multipole13[7])[i]*P13[7][m] ;	// *1
			(*Ps)[i][4][m] = (*Multipole22[10])[i]*P22[10][m] + (*Multipole22[11])[i]*P22[11][m] + (*Multipole22[12])[i]*P22[12][m] + (*Multipole13[3])[i]*P13[3][m] + (*Multipole13[4])[i]*P13[4][m];	// *b1
			(*Ps)[i][5][m] = (*Multipole22[13])[i]*P22[13][m] + (*Multipole22[14])[i]*P22[14][m] ;																										// *b2
			(*Ps)[i][6][m] = (*Multipole13[5])[i]*P13[5][m] ;																																			// *b3
			(*Ps)[i][7][m] = (*Multipole22[15])[i]*P22[15][m] + (*Multipole22[16])[i]*P22[16][m] ;																										// *b4
			(*Ps)[i][8][m] = (*Multipole22[0])[i]*P22[0][m] + (*Multipole22[3])[i]*P22[3][m] + (*Multipole22[4])[i]*P22[4][m] + (*Multipole13[0])[i]*P13[0][m] + (*Multipole13[1])[i]*P13[1][m] ;		// *b1*b1
			(*Ps)[i][9][m] = (*Multipole22[1])[i]*P22[1][m] + (*Multipole22[5])[i]*P22[5][m] ;																											// *b1*b2
			(*Ps)[i][10][m] = (*Multipole13[2])[i]*P13[2][m] ;																																			// *b1*b3
			(*Ps)[i][11][m] = (*Multipole22[2])[i]*P22[2][m] + (*Multipole22[6])[i]*P22[6][m] ;																											// *b1*b4
			(*Ps)[i][12][m] = (*Multipole22[7])[i]*P22[7][m] ;																																			// *b2*b2
			(*Ps)[i][13][m] = (*Multipole22[8])[i]*P22[8][m] ;																																			// b2*b4
			(*Ps)[i][14][m] = (*Multipole22[9])[i]*P22[9][m] ;																																			// b4*b4
		}
	}

	/***** Correlation function ******/

	// Terms c_m q^p_m in the sum of Cf11
	Eigen::ArrayXcd ArrTermsCf[Nq] ;
	for (unsigned int m = 0 ; m < Nq ; m++) {
		ArrTermsCf[m] = -ArrPow-3. ;
		ArrTermsCf[m] *= log((*s)[m]) ;
		ArrTermsCf[m] = ArrTermsCf[m].exp() ;
		ArrTermsCf[m] *= ArrCoef ;
	}

	/* Cf-linear */
	for (unsigned int i = 0 ; i < Nlmax ; i++) { // multipole
		
		Eigen::RowVectorXcd MPC11(NFFT+1) ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) 
			MPC11(m1) = MPC(2*i, -0.5*Pow[m1]) ;
		
		for (unsigned int m = 0 ; m < Nq ; m++) {
			dcomplex Cf11tmp = MPC11*ArrTermsCf[m].matrix() ;
			(*Cf)[i][0][m] = m4[i] *p.f*p.f * Cf11tmp.real() ;
			(*Cf)[i][1][m] = m2[i] *2.*p.f * Cf11tmp.real() ;
			(*Cf)[i][2][m] = m0[i] * Cf11tmp.real() ;
		}
	}

	/* Cf-counterterm */
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		
		Eigen::RowVectorXcd MPCCT(NFFT+1) ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) 
			MPCCT(m1) = MPC(2*i, -0.5*Pow[m1]-1.) ;
		
		for (unsigned int m = 0 ; m < Nq ; m++) {
			Eigen::ArrayXcd ArrTermCT ;
			ArrTermCT = -ArrPow-5. ;
			ArrTermCT *= log((*s)[m]) ;
			ArrTermCT = ArrTermCT.exp() ;
			ArrTermCT *= ArrCoef ;
			dcomplex CfCTtmp = MPCCT*ArrTermCT.matrix() ;
			(*Cf)[i][15][m] = m0[i] * 2. *CfCTtmp.real() ;			// b1*b5
			(*Cf)[i][16][m] = m2[i] * 2. *CfCTtmp.real() ;			// b1*b6
			(*Cf)[i][17][m] = m4[i] * 2. *CfCTtmp.real() ; 			// b1*b7
			(*Cf)[i][18][m] = m2[i] * 2. *p.f *CfCTtmp.real() ;		// b5
			(*Cf)[i][19][m] = m4[i] * 2. *p.f *CfCTtmp.real() ;		// b6
			(*Cf)[i][20][m] = m6[i] * 2. *p.f *CfCTtmp.real() ;		// b7
		}
	}
	
	/* Cf-1loop */
	// Copy M13 for each multipole
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		for (unsigned int j = 0 ; j < N13 ; j++)
			for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
				for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
					Mat13[i][j][m1+(NFFT+1)*m2] = Mat13[0][j][m1] ;
	}
	// Copy M22 for each multipole
	for (unsigned int i = 1 ; i < Nlmax ; i++) {
		for (unsigned int j = 0 ; j < N22 ; j++)
			for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
				for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
					Mat22[i][j][m1+(NFFT+1)*m2] = Mat22[0][j][m1+(NFFT+1)*m2] ;
	}
	// MPC x { M13, M22 }
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		// MPC from Ps_l to Cf_l
		dcomplex MPC1[NFFT+1][NFFT+1] ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
			for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
				MPC1[m1][m2] = MPC(2*i, -0.5*Pow[m1]-0.5*Pow[m2] -1.5) ;
		// MPC x M13
		for (unsigned int j = 0 ; j < N13 ; j++)
			for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
				for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
					Mat13[i][j][m1+(NFFT+1)*m2] *= MPC1[m1][m2] ;
		// MPC x M22 
		for (unsigned int j = 0 ; j < N22 ; j++)
			for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
				for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
					Mat22[i][j][m1+(NFFT+1)*m2] *= MPC1[m1][m2] ;
	}

	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		// Cf22
		double C22[N22][Nq] ;
		for (unsigned int j = 0 ; j < N22 ; j++) {
			Eigen::Map<Eigen::MatrixXcd> M22tmp(Mat22[i][j], NFFT+1, NFFT+1) ; // Eigen by default saves matrices in column-major format ; thus we used m1+(NFFT+1)*m2 indexing above
			for (unsigned int m = 0 ; m < Nq ; m++) {
				dcomplex dc22tmp = ArrTermsCf[m].matrix().transpose() * M22tmp * ArrTermsCf[m].matrix() ;
				C22[j][m] = dc22tmp.real() ;
			}
		}
		// Cf13
		double C13[N13][Nq] ;
		for (unsigned int j = 0 ; j < N13 ; j++) {
			Eigen::Map<Eigen::MatrixXcd> M13tmp(Mat13[i][j], NFFT+1, NFFT+1) ;
			for (unsigned int m = 0 ; m < Nq ; m++) {
				dcomplex dc13tmp = ArrTermsCf[m].matrix().transpose() * M13tmp * ArrTermsCf[m].matrix() ;
				C13[j][m] = dc13tmp.real() ;
			}
		}
		// Cf_l
		for (unsigned int m = 0 ; m < Nq ; m++) {
			(*Cf)[i][3][m] = (*Multipole22[17])[i]*C22[17][m] + (*Multipole22[18])[i]*C22[18][m] + (*Multipole22[19])[i]*C22[19][m] + (*Multipole13[6])[i]*C13[6][m] + (*Multipole13[7])[i]*C13[7][m] ;	// *1
			(*Cf)[i][4][m] = (*Multipole22[10])[i]*C22[10][m] + (*Multipole22[11])[i]*C22[11][m] + (*Multipole22[12])[i]*C22[12][m] + (*Multipole13[3])[i]*C13[3][m] + (*Multipole13[4])[i]*C13[4][m];	// *b1
			(*Cf)[i][5][m] = (*Multipole22[13])[i]*C22[13][m] + (*Multipole22[14])[i]*C22[14][m] ;																										// *b2
			(*Cf)[i][6][m] = (*Multipole13[5])[i]*C13[5][m] ;																																			// *b3
			(*Cf)[i][7][m] = (*Multipole22[15])[i]*C22[15][m] + (*Multipole22[16])[i]*C22[16][m] ;																										// *b4
			(*Cf)[i][8][m] = (*Multipole22[0])[i]*C22[0][m] + (*Multipole22[3])[i]*C22[3][m] + (*Multipole22[4])[i]*C22[4][m] + (*Multipole13[0])[i]*C13[0][m] + (*Multipole13[1])[i]*C13[1][m] ;		// *b1*b1
			(*Cf)[i][9][m] = (*Multipole22[1])[i]*C22[1][m] + (*Multipole22[5])[i]*C22[5][m] ;																											// *b1*b2
			(*Cf)[i][10][m] = (*Multipole13[2])[i]*C13[2][m] ;																																			// *b1*b3
			(*Cf)[i][11][m] = (*Multipole22[2])[i]*C22[2][m] + (*Multipole22[6])[i]*C22[6][m] ;																											// *b1*b4
			(*Cf)[i][12][m] = (*Multipole22[7])[i]*C22[7][m] ;																																			// *b2*b2
			(*Cf)[i][13][m] = (*Multipole22[8])[i]*C22[8][m] ;																																			// b2*b4
			(*Cf)[i][14][m] = (*Multipole22[9])[i]*C22[9][m] ;																																			// b4*b4
		}
	}

	/***** Smooth correlation function ******/

	// FFTLog decomposition
	dcomplex Coefsmooth[NFFT+1], Powsmooth[NFFT+1] ;
	CoefPow(Coefsmooth, Powsmooth, psmooth) ;

	Eigen::Map<Eigen::ArrayXcd> ArrPowsmooth(Powsmooth, NFFT+1) ;
	Eigen::Map<Eigen::ArrayXcd> ArrCoefsmooth(Coefsmooth, NFFT+1) ;

	// Terms c_m q^p_m in the sum of Cf11
	Eigen::ArrayXcd ArrTermsCfsmooth[Nq] ;
	for (unsigned int m = 0 ; m < Nq ; m++) {
		ArrTermsCfsmooth[m] = -ArrPowsmooth-3. ;
		ArrTermsCfsmooth[m] *= log((*s)[m]) ;
		ArrTermsCfsmooth[m] = ArrTermsCfsmooth[m].exp() ;
		ArrTermsCfsmooth[m] *= ArrCoefsmooth ;
	}

	/* Cf-linear */
	for (unsigned int i = 0 ; i < Nlmax ; i++) { // multipole
		
		Eigen::RowVectorXcd MPC11(NFFT+1) ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) 
			MPC11(m1) = MPC(2*i, -0.5*Powsmooth[m1]) ;
		
		for (unsigned int m = 0 ; m < Nq ; m++) {
			dcomplex Cf11tmp = MPC11*ArrTermsCfsmooth[m].matrix() ;
			(*Cfsmooth)[i][0][m] = m4[i] *p.f*p.f * Cf11tmp.real() ;
			(*Cfsmooth)[i][1][m] = m2[i] *2.*p.f * Cf11tmp.real() ;
			(*Cfsmooth)[i][2][m] = m0[i] * Cf11tmp.real() ;
		}
	}

	/* Cf-counterterm */
	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		
		Eigen::RowVectorXcd MPCCT(NFFT+1) ;
		for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) 
			MPCCT(m1) = MPC(2*i, -0.5*Pow[m1]-1.) ;
		
		for (unsigned int m = 0 ; m < Nq ; m++) {
			Eigen::ArrayXcd ArrTermCT ;
			ArrTermCT = -ArrPowsmooth-5. ;
			ArrTermCT *= log((*s)[m]) ;
			ArrTermCT = ArrTermCT.exp() ;
			ArrTermCT *= ArrCoefsmooth ;
			dcomplex CfCTtmp = MPCCT*ArrTermCT.matrix() ;
			(*Cfsmooth)[i][15][m] = m0[i] * 2. *CfCTtmp.real() ;			// b1*b5
			(*Cfsmooth)[i][16][m] = m2[i] * 2. *CfCTtmp.real() ;			// b1*b6
			(*Cfsmooth)[i][17][m] = m4[i] * 2. *CfCTtmp.real() ; 			// b1*b7
			(*Cfsmooth)[i][18][m] = m2[i] * 2. *p.f *CfCTtmp.real() ;		// b5
			(*Cfsmooth)[i][19][m] = m4[i] * 2. *p.f *CfCTtmp.real() ;		// b6
			(*Cfsmooth)[i][20][m] = m6[i] * 2. *p.f *CfCTtmp.real() ;		// b7
		}
	}

	for (unsigned int i = 0 ; i < Nlmax ; i++) {
		// Cf22
		double C22smooth[N22][Nq] ;
		for (unsigned int j = 0 ; j < N22 ; j++) {
			Eigen::Map<Eigen::MatrixXcd> M22tmp(Mat22[i][j], NFFT+1, NFFT+1) ; // Eigen by default saves matrices in column-major format ; thus we used m1+(NFFT+1)*m2 indexing above
			for (unsigned int m = 0 ; m < Nq ; m++) {
				dcomplex dc22tmp = ArrTermsCfsmooth[m].matrix().transpose() * M22tmp * ArrTermsCfsmooth[m].matrix() ;
				C22smooth[j][m] = dc22tmp.real() ;
			}
		}
		// Cf13
		double C13smooth[N13][Nq] ;
		for (unsigned int j = 0 ; j < N13 ; j++) {
			Eigen::Map<Eigen::MatrixXcd> M13tmp(Mat13[i][j], NFFT+1, NFFT+1) ;
			for (unsigned int m = 0 ; m < Nq ; m++) {
				dcomplex dc13tmp = ArrTermsCfsmooth[m].matrix().transpose() * M13tmp * ArrTermsCfsmooth[m].matrix() ;
				C13smooth[j][m] = dc13tmp.real() ;
			}
		}
		// Cf_l
		for (unsigned int m = 0 ; m < Nq ; m++) {
			(*Cfsmooth)[i][3][m] = (*Multipole22[17])[i]*C22smooth[17][m] + (*Multipole22[18])[i]*C22smooth[18][m] + (*Multipole22[19])[i]*C22smooth[19][m] + (*Multipole13[6])[i]*C13smooth[6][m] + (*Multipole13[7])[i]*C13smooth[7][m] ;	// *1
			(*Cfsmooth)[i][4][m] = (*Multipole22[10])[i]*C22smooth[10][m] + (*Multipole22[11])[i]*C22smooth[11][m] + (*Multipole22[12])[i]*C22smooth[12][m] + (*Multipole13[3])[i]*C13smooth[3][m] + (*Multipole13[4])[i]*C13smooth[4][m];	// *b1
			(*Cfsmooth)[i][5][m] = (*Multipole22[13])[i]*C22smooth[13][m] + (*Multipole22[14])[i]*C22smooth[14][m] ;																										// *b2
			(*Cfsmooth)[i][6][m] = (*Multipole13[5])[i]*C13smooth[5][m] ;																																			// *b3
			(*Cfsmooth)[i][7][m] = (*Multipole22[15])[i]*C22smooth[15][m] + (*Multipole22[16])[i]*C22smooth[16][m] ;																										// *b4
			(*Cfsmooth)[i][8][m] = (*Multipole22[0])[i]*C22smooth[0][m] + (*Multipole22[3])[i]*C22smooth[3][m] + (*Multipole22[4])[i]*C22smooth[4][m] + (*Multipole13[0])[i]*C13smooth[0][m] + (*Multipole13[1])[i]*C13smooth[1][m] ;		// *b1*b1
			(*Cfsmooth)[i][9][m] = (*Multipole22[1])[i]*C22smooth[1][m] + (*Multipole22[5])[i]*C22smooth[5][m] ;																											// *b1*b2
			(*Cfsmooth)[i][10][m] = (*Multipole13[2])[i]*C13smooth[2][m] ;																																			// *b1*b3
			(*Cfsmooth)[i][11][m] = (*Multipole22[2])[i]*C22smooth[2][m] + (*Multipole22[6])[i]*C22smooth[6][m] ;																											// *b1*b4
			(*Cfsmooth)[i][12][m] = (*Multipole22[7])[i]*C22smooth[7][m] ;																																			// *b2*b2
			(*Cfsmooth)[i][13][m] = (*Multipole22[8])[i]*C22smooth[8][m] ;																																			// b2*b4
			(*Cfsmooth)[i][14][m] = (*Multipole22[9])[i]*C22smooth[9][m] ;																																			// b4*b4
		}
	}
}