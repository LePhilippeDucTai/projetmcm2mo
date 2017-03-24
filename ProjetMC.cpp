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

int main(){

	double alpha = 0.975;
	double gamma0 =  2. ;
	int dimension = 1 ;
	double error = 1e-6;

	GaussianVector G(dimension);
	StochasticGradient S(G,gamma0,gamma,identite,alpha) ;
	S.Iterate(error);
	S.display();

	GaussianVector H(dimension);
	StochasticGradient T(H,gamma0,gamma,identite,alpha) ;
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





