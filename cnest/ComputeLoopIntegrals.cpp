#include "CBiRd.h"

#include "Eigen/Dense"

void ComputeLoopIntegrals (const ParamsP11 & p, PowerSpectraNoResum * Ps) {

        dcomplex Coef[NFFT+1], Pow[NFFT+1] ;
        CoefPow(Coef, Pow, p) ;

        Eigen::Map<Eigen::ArrayXcd> ArrPow(Pow, NFFT+1) ;
        Eigen::Map<Eigen::ArrayXcd> ArrCoef(Coef, NFFT+1) ;

        // Terms c_m k^p_m in the sum from the FFT of P11
        Eigen::ArrayXcd ArrTerms[Nk] ;
        for (unsigned int i = 0 ; i < Nk ; i++) {
                ArrTerms[i] = ArrPow ;
                ArrTerms[i] *= log(klist[i]) ;
                ArrTerms[i] = ArrTerms[i].exp() ;
                ArrTerms[i] *= ArrCoef ;
        }

        // Shot-noise 
        Eigen::ArrayXcd ArrTerm0 ;
        ArrTerm0 = ArrPow ;
        ArrTerm0 *= log(k0) ;
        ArrTerm0 = ArrTerm0.exp() ;
        ArrTerm0 *= ArrCoef ;

        // Auxiliary function involving gamma functions
        dcomplex In1n2[NFFT+1][NFFT+1] ;
        for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
                for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
                        In1n2[m1][m2] = I(-0.5*Pow[m1], -0.5*Pow[m2]) ;

        double P22[N22][Nk], P13[N13][Nk] ;

        // P22
        for (unsigned int j = 0 ; j < N22 ; j++) {

                Eigen::MatrixXcd Mat22(NFFT+1,NFFT+1) ;
                for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++)
                        for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
                                Mat22(m1,m2) = M22(j, -0.5*Pow[m1], -0.5*Pow[m2], In1n2[m1][m2], p.f) ;

                dcomplex tmp0 = ArrTerm0.matrix().transpose() * Mat22 * ArrTerm0.matrix() ;
                double ShotNoise = pow(k0,3)*tmp0.real() ;

                for (unsigned int i = 0 ; i < Nk ; i++) {
                        dcomplex tmp = ArrTerms[i].matrix().transpose() * Mat22 * ArrTerms[i].matrix() ;
                        P22[j][i] = pow(klist[i],3)*tmp.real() ;
                        P22[j][i] -= ShotNoise ;
                }
        }
        

        // P13
        for (unsigned int j = 0 ; j < N13 ; j++) {

                Eigen::RowVectorXcd Mat13(NFFT+1) ;
                for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) Mat13(m1) = M13(j, -0.5*Pow[m1], p.f) ;

                for (unsigned int i = 0 ; i < Nk ; i++) {
                        dcomplex tmp = Mat13*ArrTerms[i].matrix() ;
                        P13[j][i] = pow(klist[i],3)*P11(klist[i],p)*tmp.real() ;
                }
        }
	
	///////////////////////////////////////////
	for (unsigned int i = 0 ; i < Nl ; i++) {
    	for (unsigned int m = 0 ; m < Nk ; m++) {
    		(*Ps)[i][0][m] = (*Multipole22[17])[i]*P22[17][m] + (*Multipole22[18])[i]*P22[18][m] + (*Multipole22[19])[i]*P22[19][m] + (*Multipole13[6])[i]*P13[6][m] + (*Multipole13[7])[i]*P13[7][m] ;	// *1
    		(*Ps)[i][1][m] = (*Multipole22[10])[i]*P22[10][m] + (*Multipole22[11])[i]*P22[11][m] + (*Multipole22[12])[i]*P22[12][m] + (*Multipole13[3])[i]*P13[3][m] + (*Multipole13[4])[i]*P13[4][m];	// *b1
    		(*Ps)[i][2][m] = (*Multipole22[13])[i]*P22[13][m] + (*Multipole22[14])[i]*P22[14][m] ;																										// *b2
    		(*Ps)[i][3][m] = (*Multipole13[5])[i]*P13[5][m] ;																																			// *b3
    		(*Ps)[i][4][m] = (*Multipole22[15])[i]*P22[15][m] + (*Multipole22[16])[i]*P22[16][m] ;																										// *b4
    		(*Ps)[i][5][m] = (*Multipole22[0])[i]*P22[0][m] + (*Multipole22[3])[i]*P22[3][m] + (*Multipole22[4])[i]*P22[4][m] + (*Multipole13[0])[i]*P13[0][m] + (*Multipole13[1])[i]*P13[1][m] ;		// *b1*b1
    		(*Ps)[i][6][m] = (*Multipole22[1])[i]*P22[1][m] + (*Multipole22[5])[i]*P22[5][m] ;																											// *b1*b2
    		(*Ps)[i][7][m] = (*Multipole13[2])[i]*P13[2][m] ;																																			// *b1*b3
    		(*Ps)[i][8][m] = (*Multipole22[2])[i]*P22[2][m] + (*Multipole22[6])[i]*P22[6][m] ;																											// *b1*b4
    		(*Ps)[i][9][m] = (*Multipole22[7])[i]*P22[7][m] ;																																			// *b2*b2
    		(*Ps)[i][10][m] = (*Multipole22[8])[i]*P22[8][m] ;																																			// b2*b4
    		(*Ps)[i][11][m] = (*Multipole22[9])[i]*P22[9][m] ;																																			// b4*b4

    	}
	}

}

/* The DUMB and SLOW way

	double P22[N22][Nk], P13[N13][Nk] ;
	dcomplex P22c[N22][Nk], P13c[N22][Nk] ;

	// Terms c_m k^p_m in the sum from the FFT of P11
	dcomplex v[NFFT+1][Nk] ;
	for (unsigned int m = 0 ; m < NFFT+1 ; m++)
		for (unsigned int i = 0 ; i < Nk ; i++) 
			v[m][i] = Coef[m]*pow(klist[i],Pow[m]) ;

	// Auxiliary function involving gamma functions
	dcomplex In1n2[NFFT+1][NFFT+1] ;
	for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) 
		for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++) 
			In1n2[m1][m2] = I(-0.5*Pow[m1], -0.5*Pow[m2]) ;


	for (unsigned int m1 = 0 ; m1 < NFFT+1 ; m1++) 	{

			// P22
			for (unsigned int m2 = 0 ; m2 < NFFT+1 ; m2++)
				for (unsigned int j = 0 ; j < N22 ; j++)
					for (unsigned int i = 0 ; i < Nk ; i++)
					P22c[j][i] += v[m1][i] * M22(j, -0.5*Pow[m1], -0.5*Pow[m2], In1n2[m1][m2], p.f) * v[m2][i] ;

			// P13
			for (unsigned int j = 0 ; j < N13 ; j++)
				for (unsigned int i = 0 ; i < Nk ; i++)
					P13c[j][i] += v[m1][i] * M13(j, -0.5*Pow[m1], p.f) ;
	}


	for (unsigned int j = 0 ; j < N22 ; j++) 
		for (unsigned int i = 0 ; i < Nk ; i++)
			P22[j][i] = pow(klist[i],3)*P22c[j][i].real() ;

	for (unsigned int j = 0 ; j < N13 ; j++) 
		for (unsigned int i = 0 ; i < Nk ; i++)
			P13[j][i] = pow(klist[i],3)*P11(klist[i],p)*P13c[j][i].real() ;
	*/
