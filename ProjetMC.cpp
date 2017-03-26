//#include "ProjetMC.hpp"
#include "random_singleton.h"
#include "random_variable.hpp"
#include "random_vector.hpp"
#include "stochastic_gradient.hpp"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace arma ;

double gamma(int p,double g0){
	return g0/((double)p+100) ;
}

double identite(double x){ return x ;}

// fonction phi, identité en dimension 1
double identite(arma::vec &x){
	return x(0);
}

double sum_square(arma::vec &x){
	double res = arma::norm(x, 2);
	return res*res ;
}

double portefeuille(arma::vec &x) {
	double T = 0.25 ;
	double r = 0.05 ;
	double vol = 0.2 ;
	double K1 = 130 ;
	double K2 = 110 ;
	double S_0 = 120 ;
	double P = 0. ;
	double C = S_0*exp( (r-vol*vol*0.5)*T ) ;
	double k1 = K1/C ;
	double k2 = K2/C ;
	double S ;
	for(auto G : x){
		S = exp(vol*sqrt(T)*G) ;
		P = P + (fmax( S - k1,0.) + fmax(k2 - S,0.));
	}
	return 10.*C*P ;
}

int main(){

	double alpha = 0.995;
	double gamma0 =  1. ;
	int dimension = 5 ;
	double error = 1e-6;

	std::cout << std::endl ;
	std::cout << "Stochastic Gradient without Importance sampling : " << std::endl ;
	GaussianVector G(dimension);
	StochasticGradient S(G,gamma0,gamma,portefeuille,alpha) ;
	S.Iterate(error);
	S.display();

	std::cout << std::endl ;
	std::cout << "Stochastic Gradient with Importance sampling : " << std::endl ;
	GaussianVector H(dimension);
	StochasticGradient T(H,gamma0,gamma,portefeuille,alpha) ;
	T.IterateIS(error);
	T.display();

	return 0 ;
	//seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::cout << seed << std::endl;
	//seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::cout << seed << std::endl; 
	//Gaussian G(1.,.5);
	//std::vector<double> result = MonteCarlo(1e5,G,identite);
	//std::cout << result[0] << ' ' << result[1] << ' ' << result[2] << std::endl;

}





