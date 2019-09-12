#include "CBiRd.h"

double IntegrandBispectrumAP (const int & id, const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, const ParamsP11 & p, const double & nbar, const double & aperp, const double & apar) {
	
	double f1 = p.f ; 

	double q1, q2, q3, nu1, nu2, nu3 ;
	
	//ScoccimarroTransform (k1, k2, k3, mu, phi, q1, q2, q3, nu1, nu2, nu3) ;
	//double APfactor = 1. ;

	ScoccimarroTransformWithAP (aperp, apar, k1, k2, k3, mu, phi, q1, q2, q3, nu1, nu2, nu3) ;
	double APfactor = 1./(pow(apar,2)*pow(aperp,4)) ;

	switch (id) { 
		case 0:
			return APfactor * ( (pow(f1,3)*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(-(P11(q1,p)*pow(nu1,2)*(P11(q2,p)*pow(nu2,2)*pow(nu1*q1 + nu2*q2,2)*(3*pow(q1,4) + 3*pow(q2,4) - 14*f1*nu1*nu2*q1*q2*pow(q3,2) + pow(q2,2)*pow(q3,2) + pow(q1,2)*(-6*pow(q2,2) + pow(q3,2)) - 4*pow(q3,4)) + P11(q3,p)*pow(nu3,2)*(3*pow(q1,4) - 14*f1*nu1*nu3*q1*q3*pow(q2,2) - 4*pow(q2,4) + pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(q2,2)*pow(q3,2) + 3*pow(q3,4))*pow(nu1*q1 + nu3*q3,2))) - P11(q2,p)*P11(q3,p)*pow(nu2,2)*pow(nu3,2)*pow(nu2*q2 + nu3*q3,2)*(-4*pow(q1,4) + pow(q1,2)*(-14*f1*nu2*nu3*q2*q3 + pow(q2,2) + pow(q3,2)) + 3*pow(pow(q2,2) - pow(q3,2),2))))/14. ) ;
			break ;

		case 1:
			return APfactor * ( (pow(f1,2)*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(-(P11(q1,p)*(P11(q2,p)*(6*nu1*nu2*q2*(pow(nu1,2) + pow(nu2,2))*pow(q1,5) + 3*pow(nu1,2)*(pow(nu1,2) + pow(nu2,2))*pow(q1,6) + pow(nu2,2)*pow(q2,2)*(pow(q2,2) - pow(q3,2))*(3*(pow(nu1,2) + pow(nu2,2))*pow(q2,2) + (11*pow(nu1,2) + 4*pow(nu2,2))*pow(q3,2)) + pow(q1,4)*(-3*(2*pow(nu1,4) + pow(nu1,2)*pow(nu2,2) - pow(nu2,4))*pow(q2,2) + pow(nu1,2)*(pow(nu1,2) + 8*pow(nu2,2))*pow(q3,2)) - 2*nu1*nu2*q2*pow(q1,3)*(6*(pow(nu1,2) + pow(nu2,2))*pow(q2,2) + (7*f1*pow(nu1,4) - pow(nu2,2) + pow(nu1,2)*(-1 + 14*f1*pow(nu2,2)))*pow(q3,2)) + 2*nu1*nu2*q1*q2*(3*(pow(nu1,2) + pow(nu2,2))*pow(q2,4) + (pow(nu1,2) + pow(nu2,2) - 14*f1*pow(nu1,2)*pow(nu2,2) - 7*f1*pow(nu2,4))*pow(q2,2)*pow(q3,2) - 4*(pow(nu1,2) + pow(nu2,2))*pow(q3,4)) + pow(q1,2)*(3*(pow(nu1,4) - pow(nu1,2)*pow(nu2,2) - 2*pow(nu2,4))*pow(q2,4) + (pow(nu1,4)*(1 - 42*f1*pow(nu2,2)) + pow(nu2,4) + pow(nu1,2)*(16*pow(nu2,2) - 42*f1*pow(nu2,4)))*pow(q2,2)*pow(q3,2) - pow(nu1,2)*(4*pow(nu1,2) + 11*pow(nu2,2))*pow(q3,4))) + P11(q3,p)*(6*nu1*nu3*q3*(pow(nu1,2) + pow(nu3,2))*pow(q1,5) + 3*pow(nu1,2)*(pow(nu1,2) + pow(nu3,2))*pow(q1,6) + pow(nu3,2)*pow(q3,2)*(-pow(q2,2) + pow(q3,2))*((11*pow(nu1,2) + 4*pow(nu3,2))*pow(q2,2) + 3*(pow(nu1,2) + pow(nu3,2))*pow(q3,2)) - 2*nu1*nu3*q3*pow(q1,3)*((7*f1*pow(nu1,4) - pow(nu3,2) + pow(nu1,2)*(-1 + 14*f1*pow(nu3,2)))*pow(q2,2) + 6*(pow(nu1,2) + pow(nu3,2))*pow(q3,2)) + pow(q1,4)*((pow(nu1,4) + 8*pow(nu1,2)*pow(nu3,2))*pow(q2,2) + 3*(-2*pow(nu1,4) - pow(nu1,2)*pow(nu3,2) + pow(nu3,4))*pow(q3,2)) + 2*nu1*nu3*q1*q3*(-4*(pow(nu1,2) + pow(nu3,2))*pow(q2,4) + (pow(nu1,2) + pow(nu3,2) - 14*f1*pow(nu1,2)*pow(nu3,2) - 7*f1*pow(nu3,4))*pow(q2,2)*pow(q3,2) + 3*(pow(nu1,2) + pow(nu3,2))*pow(q3,4)) + pow(q1,2)*(-((4*pow(nu1,4) + 11*pow(nu1,2)*pow(nu3,2))*pow(q2,4)) + (pow(nu1,4)*(1 - 42*f1*pow(nu3,2)) + pow(nu3,4) + pow(nu1,2)*(16*pow(nu3,2) - 42*f1*pow(nu3,4)))*pow(q2,2)*pow(q3,2) + 3*(pow(nu1,4) - pow(nu1,2)*pow(nu3,2) - 2*pow(nu3,4))*pow(q3,4))))) + P11(q2,p)*P11(q3,p)*(pow(q1,4)*(8*nu2*nu3*q2*q3*(pow(nu2,2) + pow(nu3,2)) + (4*pow(nu2,4) + 11*pow(nu2,2)*pow(nu3,2))*pow(q2,2) + pow(nu3,2)*(11*pow(nu2,2) + 4*pow(nu3,2))*pow(q3,2)) - pow(q1,2)*(2*nu2*nu3*q3*(-7*f1*pow(nu2,4) + pow(nu3,2) + pow(nu2,2)*(1 - 14*f1*pow(nu3,2)))*pow(q2,3) + (pow(nu2,4) + 8*pow(nu2,2)*pow(nu3,2))*pow(q2,4) + (pow(nu2,4)*(1 - 42*f1*pow(nu3,2)) + pow(nu3,4) + pow(nu2,2)*(16*pow(nu3,2) - 42*f1*pow(nu3,4)))*pow(q2,2)*pow(q3,2) + 2*nu2*nu3*q2*(pow(nu3,2) + pow(nu2,2)*(1 - 14*f1*pow(nu3,2)) - 7*f1*pow(nu3,4))*pow(q3,3) + pow(nu3,2)*(8*pow(nu2,2) + pow(nu3,2))*pow(q3,4)) - 3*(pow(nu2,2) + pow(nu3,2))*pow(nu2*q2 + nu3*q3,2)*pow(pow(q2,2) - pow(q3,2),2))))/14. ) ;
			break ;

		case 2:
			return APfactor * ( (pow(f1,2)*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q2,p)*P11(q3,p)*pow(nu2,2)*pow(nu3,2)*pow(q1,2)*(pow(q1,4) + pow(q2,4) + 12*pow(q2,2)*pow(q3,2) - 2*pow(q1,2)*(pow(q2,2) + pow(q3,2)) + pow(q3,4)) + P11(q1,p)*pow(nu1,2)*(P11(q3,p)*pow(nu3,2)*pow(q2,2)*(pow(q1,4) - 2*pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)) + P11(q2,p)*pow(nu2,2)*pow(q3,2)*(pow(q1,4) + 2*pow(q1,2)*(6*pow(q2,2) - pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)))))/7. ) ;
			break ;

		case 3:
			return APfactor * ( 2*pow(f1,2)*(P11(q2,p)*P11(q3,p)*pow(nu2,2)*pow(nu3,2) + P11(q1,p)*pow(nu1,2)*(P11(q2,p)*pow(nu2,2) + P11(q3,p)*pow(nu3,2))) ) ;
			break ;

		case 4:
			return APfactor * ( (P11(q1,p) + P11(q2,p) + P11(q3,p))*pow(nbar,-1) ) ;
			break ;

		case 5:
			return APfactor * ( (f1*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q1,p)*(P11(q2,p)*(-6*nu1*nu2*q2*pow(q1,5) - 3*pow(nu1,2)*pow(q1,6) + pow(q1,4)*((6*pow(nu1,2) - 3*pow(nu2,2))*pow(q2,2) - (8*pow(nu1,2) + 7*pow(nu2,2))*pow(q3,2)) - pow(q2,2)*(pow(q2,2) - pow(q3,2))*(3*pow(nu2,2)*pow(q2,2) + (7*pow(nu1,2) + 11*pow(nu2,2))*pow(q3,2)) + 2*nu1*nu2*q2*pow(q1,3)*(6*pow(q2,2) + (-1 + 7*f1*(2*pow(nu1,2) + pow(nu2,2)))*pow(q3,2)) + 2*nu1*nu2*q1*q2*(-3*pow(q2,4) + (-1 + 7*f1*(pow(nu1,2) + 2*pow(nu2,2)))*pow(q2,2)*pow(q3,2) + 4*pow(q3,4)) + pow(q1,2)*(-3*(pow(nu1,2) - 2*pow(nu2,2))*pow(q2,4) + (14*f1*pow(nu1,4) + pow(nu2,2)*(-15 + 14*f1*pow(nu2,2)) + pow(nu1,2)*(-15 + 56*f1*pow(nu2,2)))*pow(q2,2)*pow(q3,2) + (11*pow(nu1,2) + 7*pow(nu2,2))*pow(q3,4))) + P11(q3,p)*(-6*nu1*nu3*q3*pow(q1,5) - 3*pow(nu1,2)*pow(q1,6) + 2*nu1*nu3*q3*pow(q1,3)*((-1 + 7*f1*(2*pow(nu1,2) + pow(nu3,2)))*pow(q2,2) + 6*pow(q3,2)) + (pow(q2,2) - pow(q3,2))*pow(q3,2)*((7*pow(nu1,2) + 11*pow(nu3,2))*pow(q2,2) + 3*pow(nu3,2)*pow(q3,2)) - pow(q1,4)*((8*pow(nu1,2) + 7*pow(nu3,2))*pow(q2,2) + 3*(-2*pow(nu1,2) + pow(nu3,2))*pow(q3,2)) + 2*nu1*nu3*q1*q3*(4*pow(q2,4) + (-1 + 7*f1*(pow(nu1,2) + 2*pow(nu3,2)))*pow(q2,2)*pow(q3,2) - 3*pow(q3,4)) + pow(q1,2)*((11*pow(nu1,2) + 7*pow(nu3,2))*pow(q2,4) + (14*f1*pow(nu1,4) + pow(nu3,2)*(-15 + 14*f1*pow(nu3,2)) + pow(nu1,2)*(-15 + 56*f1*pow(nu3,2)))*pow(q2,2)*pow(q3,2) - 3*(pow(nu1,2) - 2*pow(nu3,2))*pow(q3,4)))) + P11(q2,p)*P11(q3,p)*(pow(q1,4)*(8*nu2*nu3*q2*q3 + (11*pow(nu2,2) + 7*pow(nu3,2))*pow(q2,2) + (7*pow(nu2,2) + 11*pow(nu3,2))*pow(q3,2)) + pow(q1,2)*(2*nu2*nu3*q3*(-1 + 7*f1*(2*pow(nu2,2) + pow(nu3,2)))*pow(q2,3) - (8*pow(nu2,2) + 7*pow(nu3,2))*pow(q2,4) + (14*f1*pow(nu2,4) + pow(nu3,2)*(-15 + 14*f1*pow(nu3,2)) + pow(nu2,2)*(-15 + 56*f1*pow(nu3,2)))*pow(q2,2)*pow(q3,2) + 2*nu2*nu3*q2*(-1 + 7*f1*(pow(nu2,2) + 2*pow(nu3,2)))*pow(q3,3) - (7*pow(nu2,2) + 8*pow(nu3,2))*pow(q3,4)) - 3*pow(nu2*q2 + nu3*q3,2)*pow(pow(q2,2) - pow(q3,2),2))))/14. ) ;
			break ;

		case 6:
			return APfactor * ( (f1*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q2,p)*P11(q3,p)*(pow(nu2,2) + pow(nu3,2))*pow(q1,2)*(pow(q1,4) + pow(q2,4) + 12*pow(q2,2)*pow(q3,2) - 2*pow(q1,2)*(pow(q2,2) + pow(q3,2)) + pow(q3,4)) + P11(q1,p)*(P11(q3,p)*(pow(nu1,2) + pow(nu3,2))*pow(q2,2)*(pow(q1,4) - 2*pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)) + P11(q2,p)*(pow(nu1,2) + pow(nu2,2))*pow(q3,2)*(pow(q1,4) + 2*pow(q1,2)*(6*pow(q2,2) - pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)))))/7. ) ;
			break ;

		case 7:
			return APfactor * ( 2*f1*(P11(q2,p)*P11(q3,p)*(pow(nu2,2) + pow(nu3,2)) + P11(q1,p)*(P11(q2,p)*(pow(nu1,2) + pow(nu2,2)) + P11(q3,p)*(pow(nu1,2) + pow(nu3,2)))) ) ;
			break ;

		case 8:
			return APfactor * ( (pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q1,p)*(P11(q2,p)*pow(q3,2)*(2*f1*nu1*nu2*q2*pow(q1,3) - pow(q1,4) + 2*f1*nu1*nu2*q1*pow(q2,3) - pow(q2,4) + pow(q2,2)*pow(q3,2) + pow(q1,2)*(2*(-1 + f1*(pow(nu1,2) + pow(nu2,2)))*pow(q2,2) + pow(q3,2))) + P11(q3,p)*pow(q2,2)*(2*f1*nu1*nu3*q3*pow(q1,3) - pow(q1,4) + (pow(q2,2) - pow(q3,2))*pow(q3,2) + pow(q1,2)*(pow(q2,2) + 2*(-1 + f1*(pow(nu1,2) + pow(nu3,2)))*pow(q3,2)) + 2*f1*nu1*nu3*q1*pow(q3,3))) + P11(q2,p)*P11(q3,p)*pow(q1,2)*(2*f1*nu2*nu3*q3*pow(q2,3) - pow(q2,4) + 2*(-1 + f1*(pow(nu2,2) + pow(nu3,2)))*pow(q2,2)*pow(q3,2) + pow(q1,2)*(pow(q2,2) + pow(q3,2)) + 2*f1*nu2*nu3*q2*pow(q3,3) - pow(q3,4))))/2. ) ;
			break ;

		case 9:
			return APfactor * ( (pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q2,p)*P11(q3,p)*pow(q1,2)*(pow(q1,4) + pow(q2,4) + 12*pow(q2,2)*pow(q3,2) - 2*pow(q1,2)*(pow(q2,2) + pow(q3,2)) + pow(q3,4)) + P11(q1,p)*(P11(q3,p)*pow(q2,2)*(pow(q1,4) - 2*pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)) + P11(q2,p)*pow(q3,2)*(pow(q1,4) + 2*pow(q1,2)*(6*pow(q2,2) - pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)))))/7. ) ;
			break ;

		case 10:
			return APfactor * ( 2*(P11(q2,p)*P11(q3,p) + P11(q1,p)*(P11(q2,p) + P11(q3,p))) ) ;
			break ;

	}


}
void ScoccimarroTransform (const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) {

	double cos12 = -(k1*k1 + k2*k2 - k3*k3) / (2.*k1*k2) ;

	nu1 = mu ;
	nu2 = mu * cos12 - sqrt(1-mu*mu)*sqrt(1.-cos12*cos12)*cos(phi) ;
	nu3 = - k1/k3 * nu1 - k2/k3 * nu2 ;

	q1 = k1 ;
	q2 = k2 ;
	q3 = k3 ;
}


// If the condition of closed triangle is broken by the AP effect, return 0 and the integrands are set to 0.
void ScoccimarroTransformWithAP (const double & aperp, const double & apar, 
	const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) {

	double F = apar/aperp ;

	double cos12 = -(k1*k1 + k2*k2 - k3*k3) / (2.*k1*k2) ;

	nu1 = mu ;
	nu2 = mu * cos12 - sqrt(1-mu*mu)*sqrt(1.-cos12*cos12)*cos(phi) ;
	nu3 = - k1/k3 * nu1 - k2/k3 * nu2 ;

	q1 = k1/aperp * sqrt(1.+ nu1*nu1 * (pow(F,-2)-1.)) ;
	q2 = k2/aperp * sqrt(1.+ nu2*nu2 * (pow(F,-2)-1.)) ;
	q3 = k3/aperp * sqrt(1.+ nu3*nu3 * (pow(F,-2)-1.)) ;

	nu1 = nu1/F * 1./sqrt(1.+ nu1*nu1 * (pow(F,-2)-1.)) ;
	nu2 = nu2/F * 1./sqrt(1.+ nu2*nu2 * (pow(F,-2)-1.)) ;
	nu3 = nu3/F * 1./sqrt(1.+ nu3*nu3 * (pow(F,-2)-1.)) ;
	
}