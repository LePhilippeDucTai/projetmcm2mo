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

template <typename T> 
void display(std::vector<T> &V){
	for (auto value : V) {
    std::cout << value << ' ' ;
	}
	std::cout<< std::endl ;;
}

int main(){

	Random::Randomize(100);

	double alpha = 0.975 ;
	double xi0 = 0. ;
	double gamma0 =  2. ;

	StochasticGradient S(gamma0,gamma,xi0,alpha) ;
	double error = 1e-7;

	int i = 1 ;
	double N ;

	mat V = randu<mat>(5,6);
	V.print("test: ");
	//seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::cout << seed << std::endl;
	//seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::cout << seed << std::endl; 
	//Gaussian G(1.,.5);
	//std::vector<double> result = MonteCarlo(1e5,G,identite);
	//std::cout << result[0] << ' ' << result[1] << ' ' << result[2] << std::endl;
	
	UniformVector U(0,1,10);
	U.init(100);
	mat X = U();
	X.print("Uniform Vector :");

	int d = 5;
	 vec mu = zeros<vec>(d);
	 mat sigma = eye<mat>(d,d);
	mu.print("test:");
	sigma.print("test2:");

	// GaussianVector G(mu,sigma,10);
	 GaussianVector G(3);
	G.init(10);
	mat C = G() ;
	C.print("Gaussian Vector :");


 	vec mu1 = {1,2,3};
 	mat mu2 = conv_to<mat>::from(mu1);
 	mu2.print("mu2: ");
 	mat sigma1 = {  {1.,0.9, 0.1},
 					{0.9, 1., 0.5},
 					{0.1, 0.5, 1.1}};
 	mat A = chol(sigma1);
 	A.print("Cholesky decomp:");


 	GaussianVector G1(mu1,sigma1,3);
	G1.init(10);
	mat C1 = G1() ;
	C1.print("Gaussian Vector :");
 	mu1.print("Mu : ");
 	sigma1.print("Sigma :");

 	std::cout << "Gaussian density at mu " << G1.Pdf(mu1) << std::endl;
	std::cout << "Gaussian density at mu+0.1 " << G1.Pdf(mu1+0.1) << std::endl;
	std::cout << "Gaussian density at mu-0.1 " << G1.Pdf(mu1-0.1) << std::endl;
	std::cout << "Gaussian density at mu+4 " << G1.Pdf(mu1+4) << std::endl;
	std::cout << "Gaussian density at mu-4 " << G1.Pdf(mu1-4) << std::endl;


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





