//#include "ProjetMC.hpp"
#include "random_singleton.h"
#include "stochastic_gradient.hpp"
#include <iostream>
#include <cmath>

double gamma(int p,double g0){
	return g0/((double)p+100) ;
}


int main(){

Random::Randomize(100);

double alpha = 0.95 ;
double xi0 = 0. ;
double gamma0 =  2. ;

StochasticGradient S(gamma0,gamma,xi0,alpha) ;
double error = 1e-7;

int i = 1 ;
double N ;
do {
	N = Random::Gaussian();
	//N = Random::Exponential(1);
	S.Iterate(i,N); 
	i = i+1 ;
}while( S.precision_xi() + S.precision_c() > error);
S.display();
//std::cout << (-1)*log(1-alpha) << std::endl ;
return 0 ;
}





