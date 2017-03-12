#include <iostream>
#include <random>
#include <cstdio>
#include <ctime>
#include "utility.hpp"

std::mt19937 mt_rand(time(0));

int rpois(double lambda){
	std::poisson_distribution<int> pois(lambda);
	return pois(mt_rand) ;
}

double rlnorm(double mu,double sigma){
	std::lognormal_distribution<double> ln(mu,sigma);
	return ln(mt_rand) ;
}


double PoissonComp(double lambda, double mu,double sigma){
	double X = 0 ;
	int N = rpois(lambda);
	for(int i=0 ; i < N ; i++) {
		X = X + rlnorm(mu,sigma);
	}
	return X ;
}

std::vector<double> vPoissonComp(int nsim, double lambda, double mu, double sigma) {
	std::vector<double> X ;
	for(int i=0 ; i < nsim ; i++ ) {
		X.push_back(PoissonComp(lambda,mu,sigma));
	}
	return X ;
}

int main(){
	// double lambda = 100 ;
	// double mu = 1 ;
	// double sigma = 1 ;
	// int nsim = 100000 ;

 //    std::clock_t start;
 //    double duration;

 //    start = std::clock();
	// std::vector<double> v = vPoissonComp(nsim,lambda,mu,sigma);

 //  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

 //    std::cout<<"Duration : "<< duration <<'\n';

//std::vector<double> r = range(0, 1, 0.01) ;
std::vector<double> r ;
display_vect(r);




}
