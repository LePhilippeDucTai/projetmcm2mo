#ifndef STOCHASTIC_GRADIENT_HPP
#define STOCHASTIC_GRADIENT_HPP
#include "random_singleton.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#include "random_vector.hpp"

// template <typename randomvector> 
class StochasticGradient{

public:
	// L'algorithme du gradient stochastique est déterminé par la suite gamma et son terme initial.
	//Constructeur
	// StochasticGradient(double gamma0, double (*gamma)(int,double), double xi0,double alpha);
	StochasticGradient(GaussianVector X, double gamma0, double (*gamma)(int,double), double (*phi)(const arma::vec &),double alpha);
	
	// Initialize the sequences before doing iterations
	void init_seq(double xi0, double c0, arma::vec theta, arma::vec mu,double rho, double a, double b);

	// Calcul des variables _xi et _c par itération
	//double gamma(int n);
	void Iterate(double epsilon); // L'algorithme s'arrête quand la convergence est suffisamment fine
	void IterateIS(double epsilon);
	inline double precision_xi() const;
	inline double precision_c() const ;
	double H1(const arma::vec &x);
	double H2(const arma::vec &x);
	double L1(const arma::vec &x);
	double L2(const arma::vec &x);
	arma::vec K1(arma::vec &x);
	arma::vec K2(arma::vec &x);

	// Accesseurs
	double getVaR();
	double getCVaR();
	void display();
private:
	double _xi, _lastxi ; // représente la VaR
	double _c, _lastc  ; // représente la CVaR
	double (*_gamma)(int,double) ; //Fonction gamma
	double _gamma0 ; // Terme initial de la suite gamma
	double (*_phi)(const arma::vec &x); // Fonction "payoff" de R^d -> R
	double _alpha ;
	int niterate ;
	arma::vec _theta ;
	arma::vec _mu ; 
	double _b ; // = 2 in gaussian case
	double _rho ; // = 1 in gaussian case
	double _a ;
	int _dim ; // Dimension de X

	GaussianVector X ; // Find a way to generalize to RandomVectors

};

#endif
