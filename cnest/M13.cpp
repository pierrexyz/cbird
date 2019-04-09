#include "CBiRd.h"

dcomplex M13 (const unsigned int & j, const dcomplex & n1, const double & f1) {
	switch(j) {
		case 0:
			return (9.*f1*pow(n1,-1)*pow(M_PI,-1)*pow(-6. + 11.*n1 - 6.*pow(n1,2) + pow(n1,3),-1)*tan(n1*M_PI))/56. ;
			break ; 
		case 1:
			return (9.*pow(n1,-1)*pow(M_PI,-1)*pow(-6. + 11.*n1 - 6.*pow(n1,2) + pow(n1,3),-1)*tan(n1*M_PI))/112. ;
			break ; 
		case 2:
			return -(pow(n1,-1)*pow(M_PI,-1)*pow(-6. + 5.*n1 + 5.*pow(n1,2) - 5.*pow(n1,3) + pow(n1,4),-1)*tan(n1*M_PI))/14. ;
			break ; 
		case 3:
			return (-3.*f1*(1. + 3.*f1 - 3.*n1)*pow(-3. + n1,-1)*pow(-2. + n1,-1)*pow(-1. + n1,-1)*pow(n1,-1)*pow(1. + n1,-1)*pow(M_PI,-1)*tan(n1*M_PI))/56. ;
			break ; 
		case 4:
			return (9.*(1. + 2.*n1)*pow(f1,2)*pow(n1,-1)*pow(M_PI,-1)*pow(-6. + 5.*n1 + 5.*pow(n1,2) - 5.*pow(n1,3) + pow(n1,4),-1)*tan(n1*M_PI))/56. ;
			break ; 
		case 5:
			return -(f1*pow(n1,-1)*pow(M_PI,-1)*pow(-6. + 5.*n1 + 5.*pow(n1,2) - 5.*pow(n1,3) + pow(n1,4),-1)*tan(n1*M_PI))/14. ;
			break ; 
		case 6:
			return (-3.*(5. + 6.*f1 - 3.*n1)*pow(f1,2)*pow(-3. + n1,-1)*pow(-2. + n1,-1)*pow(-1. + n1,-1)*pow(n1,-1)*pow(1. + n1,-1)*pow(M_PI,-1)*tan(n1*M_PI))/112. ;
			break ; 
		case 7:
			return (9.*pow(f1,3)*pow(M_PI,-1)*pow(-6. + 5.*n1 + 5.*pow(n1,2) - 5.*pow(n1,3) + pow(n1,4),-1)*tan(n1*M_PI))/56. ;
			break ; 
	}
}
