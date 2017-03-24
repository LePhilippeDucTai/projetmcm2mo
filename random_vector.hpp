#ifndef RANDOM_VECTOR_HPP
#define RANDOM_VECTOR_HPP
#include <functional>
#include <vector>
#include <cmath>
#include <list>
#include <iostream>
#include <random>
#include <chrono>
#include "random_variable.hpp"
#include <armadillo>

template <class T> struct RandomVect {
	
	typedef T result_type;
	RandomVect() 
			: generator(1) {} // default seed = 1 
	RandomVect(T value)
			: value(value),generator(1) {}

	virtual ~RandomVect() {}
	virtual T operator()() = 0; 
	virtual double Pdf(T) = 0 ;
	virtual arma::vec Gradient(T) = 0 ;
	int Size() const {
		return d;
	}
	T Current() const {
		return value;
	}
	// Initializes the seed and the dimension of the vector
	void init(int seed_) {
		seed = seed_ ;
		generator.seed(seed);
	}

	double SGIS_rho(){
		return _SGIS_rho ;
	}
	double SGIS_a(){
		return _SGIS_a ;
	}
	double SGIS_b(){
		return _SGIS_b;
	}
	double SGIS_C(){
		return _SGIS_C;
	}
 protected:
	T value;
	int seed ;
	int d ; // Dimension of the vector
	std::mt19937 generator;

	 // Specialized constants in the stochastic gradient+importance sampling algorithm
	double _SGIS_rho, _SGIS_a, _SGIS_b, _SGIS_C ;
};


struct UniformVector : public RandomVect<arma::vec> {
	UniformVector(double left_boundary , double right_boundary, int N)
			: left_boundary(left_boundary), right_boundary(right_boundary) 
			{ d = N; }
	arma::vec operator()() {
		arma::vec A(d);
		value = A ;
		value.imbue( [&]() { return left_boundary + (right_boundary - left_boundary) * 
								generator()/(static_cast<double>(generator.max()));  } );
		return value ;
	}

	inline double Pdf(arma::vec x) {
		for(int i = 0 ; i < d ; i++ ) {
			if (x(i) < left_boundary && x(i) > right_boundary)
				return 0. ;
		}
		return 1./pow(static_cast<double>(right_boundary - left_boundary),d) ;
	}
 private:
	double left_boundary, right_boundary;
};


struct GaussianVector : public RandomVect<arma::vec> {
	// Gaussian vector initialized as N(0, Id)
	 GaussianVector(int N) // flag == true if sigma == Id 
 			: mu(arma::zeros<arma::vec>(N)), sigma(arma::eye<arma::mat>(N,N)), flag_cr(true) 
				{ d = N; 
				_SGIS_rho = 1. ;
				_SGIS_a = 1. ;
				_SGIS_b = 2. ;
				_SGIS_C = 2. ;
				sigmainv = sigma ;
				}
 	GaussianVector(arma::vec mu , arma::mat sigma, int N)
			: mu(mu), sigma(sigma) {
				d = N;
	 			sigmachol = arma::chol(sigma,"lower");
	 			sigmainv = arma::inv(sigma) ; 
	 			sigmadetsqrt = sqrt(arma::det(sigma));
	 			flag_cr = false ; // flag == false if sigma != Id
	 			_SGIS_rho = 1. ;
				_SGIS_a = 1. ;
				_SGIS_b = 2. ;
				_SGIS_C = 2. ;
	 		}

	arma::vec operator()(){
		value.zeros(d);
		value.imbue([&](){ return norm_inv(generator()/(static_cast<double> (generator.max())));});
		if( flag_cr )
			return value ;
		else return value = mu + sigmachol*value ;
	}

	inline double Pdf(arma::vec x) { 
		double C = 1./(pow(2.*M_PI,d/2.)*sigmadetsqrt);
		double Y = arma::as_scalar(-0.5*(x-mu).t()*sigmainv*(x-mu)) ;
		return C*exp(Y) ;
	}

	// Gradient de la Pdf
	inline arma::vec Gradient(arma::vec x){
		return Pdf(x)*sigmainv*(x-mu) ;
	}

	inline double norm_inv(double x){
		int signe = 0;
		double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328, d1 = 1.432788, d2 = 0.189269, d3 = 0.001308 ;

		if(x>0.5){
			signe = 1 ;
			x=1.0-x;
		}
		else 
			signe=-1.0;
		double t = sqrt(-2.0*log(x));
		return (signe * (t-((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t)))));
	}

 private:
 	arma::vec mu ;
 	arma::mat sigma ;
 	arma::mat sigmainv ;
 	arma::mat sigmachol ;
 	double sigmadetsqrt ;
 	bool flag_cr ;
};











#endif
