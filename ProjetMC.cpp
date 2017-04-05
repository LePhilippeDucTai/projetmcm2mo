//#include "ProjetMC.hpp"
#include "random_singleton.h"
#include "random_variable.hpp"
#include "random_vector.hpp"
#include "stochastic_gradient.hpp"
#include <iostream>
#include <cmath>
#include <armadillo>
// #include "RNGSobol.h"

using namespace arma ;

double gamma(int p,double g0){
	return g0/((double)p+100) ;
}

double identite(double x){ return x ;}

// fonction phi, identité en dimension 1
double identite(const arma::vec &x){
	return x(0);
}

double sum_square(const arma::vec &x){
	double res = arma::norm(x, 2);
	return res*res ;
}

double portefeuille(const arma::vec &x) {
	double T = 0.25 ;
	double r = 0.05 ;
	double vol = 0.2 ;
	double K1 = 130 ;
	double K2 = 110 ;
	double S_0 = 120 ;
	double C = S_0*exp( (r-vol*vol*0.5)*T ) ;
	double k1 = K1/C ;
	double k2 = K2/C ;
	double S ;
	double P = 0. ;
	for(auto G : x){
		S = exp(vol*sqrt(T)*G) ;
		P = P + (fmax( S - k1,0.) + fmax(k2 - S,0.));
	}
	return 10.*C*P ;
}

int main(){

	double alpha = 0.975;
	double gamma0 =  1. ;
	int dimension = 1 ;
	double error = 1e-6;
	// RNGSobol X ;
	// X.setDim(1);
	// X.setSeed(1);
	// std::cout << X.generate01() << std::endl;
	std::cout << "##########################################" << std::endl ;
	alpha = 0.975;
	dimension = 1 ;
	std::cout << std::endl ;
	std::cout << "####### Application à des Gaussiennes Centrées réduites 97.5% ########" << std::endl ;
	std::cout << "Stochastic Gradient without Importance sampling : " << std::endl ;
	GaussianVector G1(dimension);
	StochasticGradient S1(G1,gamma0,gamma,identite,alpha) ;
	S1.Iterate(error);
	S1.display();

	std::cout << std::endl ;
	std::cout << "Stochastic Gradient with Importance sampling : " << std::endl ;
	GaussianVector H1(dimension);
	StochasticGradient T1(H1,gamma0,gamma,identite,alpha) ;
	T1.IterateIS(error);
	T1.display();

	std::cout << "##########################################" << std::endl ;
	alpha = 0.975;
	dimension = 5 ;
	std::cout << std::endl ;
	std::cout << "####### Application à des Chi-deux à 5 degres de liberte 97.5% ########" << std::endl ;
	std::cout << "Stochastic Gradient without Importance sampling : " << std::endl ;
	GaussianVector G2(dimension);
	StochasticGradient S2(G2,gamma0,gamma,sum_square,alpha) ;
	S2.Iterate(error);
	S2.display();

	std::cout << std::endl ;
	std::cout << "Stochastic Gradient with Importance sampling : " << std::endl ;
	GaussianVector H2(dimension);
	StochasticGradient T2(H2,gamma0,gamma,sum_square,alpha) ;
	T2.IterateIS(error);
	T2.display();




	std::cout << "##########################################" << std::endl ;

	alpha = 0.995;
	dimension = 5 ;
		std::cout << std::endl ;
	std::cout << "####### Exemple 2 : Calls et Puts 99.5% ########" << std::endl ;
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





