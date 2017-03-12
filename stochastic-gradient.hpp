#include "random-singleton.h"
#include <iostream>
#include <cmath>
#include <vector>

class StochasticGradient{

public:
	// L'algorithme du gradient stochastique est déterminé par la suite gamma et son terme initial.
	//Constructeur
	StochasticGradient(double gamma0, double (*gamma)(int,double), double xi0,double alpha);
	
	// Calcul des variables _xi et _c par itération
	//double gamma(int n);
	void Iterate(int n, double X); // Prend en entrée un entier n et X une variable aléatoire
	double precision_xi();
	double precision_c();
	double H1(double x);
	double H2(double x);
	// Accesseurs
	double getVaR();
	double getCVaR();
	void display();
private:
	double _xi, _lastxi ; // représente la VaR
	double _c, _lastc  ; // représente la CVaR
	double (*_gamma)(int,double) ; //Fonction gamma
	double _gamma0 ; // Terme initial de la suite gamma
	double _alpha ;
	int niterate ;

};