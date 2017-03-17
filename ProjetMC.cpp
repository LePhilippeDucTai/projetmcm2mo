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

double identite(arma::vec x){
	return x(0);
}

// template <typename T> 
// void display(std::vector<T> &V){
// 	for (auto value : V) {
//     std::cout << value << ' ' ;
// 	}
// 	std::cout<< std::endl ;;
// }

int main(){

	double alpha = 0.975 ;
	double gamma0 =  2. ;
	int dim = 1 ;
	double error = 1e-6;

	GaussianVector G(dim);
	G.init(10);
	mat C = G() ;
	C.print("Gaussian Vector :");

	StochasticGradient S(G,gamma0,gamma,identite,alpha,dim) ;
	S.Iterate(error);
	S.display();



	//seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::cout << seed << std::endl;
	//seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::cout << seed << std::endl; 
	//Gaussian G(1.,.5);
	//std::vector<double> result = MonteCarlo(1e5,G,identite);
	//std::cout << result[0] << ' ' << result[1] << ' ' << result[2] << std::endl;
	


 	// vec mu1 = {1,2,3};
 	// mat mu2 = conv_to<mat>::from(mu1);
 	// mu2.print("mu2: ");
 	// mat sigma1 = {  {1.,0.9, 0.1},
 	// 				{0.9, 1., 0.5},
 	// 				{0.1, 0.5, 1.1}};
 	// mat A = chol(sigma1);
 	// A.print("Cholesky decomp:");


 	// GaussianVector G1(mu1,sigma1,3);



	// exit(-1);
	// do {
	// 	N = Random::Gaussian();
	// 	//N = Random::Exponential(1);
	// 	S.Iterate(i,N); 
	// 	i = i+1 ;
	// } while( S.precision_xi() + S.precision_c() > error);
	// S.display();
	// //std::cout << (-1)*log(1-alpha) << std::endl ;
	// return 0 ;
}





