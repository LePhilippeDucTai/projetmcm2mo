#include "stochastic_gradient.hpp"

StochasticGradient::StochasticGradient(double gamma0,double (*gamma)(int,double), double xi0, double alpha) 
		: _gamma0(gamma0), _gamma(gamma), _xi(xi0), _alpha(alpha), _c(0.), niterate(0.) {}

// double StochasticGradient::gamma(int n){
// 	return (_gamma0)/((double)n) ;
// }

// Cette fonction sera appelée dans une boucle pour i
// Elle calculera le taux de variation de la dernière valeur de xi avec la nouvelle calculée
// On s'arrêtera dans la boucle lorsque le taux d'erreur sera inférieur à une certaine valeur
void StochasticGradient::Iterate(int n, double fX){
	_lastxi = _xi ; // Enregistrement de la dernière valeur de xi
	_lastc = _c ; // Enregistrement de la dernière valeur de c
	double h1 = StochasticGradient::H1(fX);
	double h2 = StochasticGradient::H2(fX); 
	_xi = _xi - _gamma(n,_gamma0)*h1 ;
	_c = _c - _gamma(n,_gamma0)*h2 ;
	niterate++;
}

//Calcul de la différence |xi_(n+1)-xi_(n)| ;
double StochasticGradient::precision_xi(){
	return abs(_xi - _lastxi) ;
}

double StochasticGradient::precision_c(){
	return abs(_c - _lastc) ;
}

// Définition des accesseurs
double StochasticGradient::getVaR(){
	return _xi ;
}

double StochasticGradient::getCVaR(){
	return _c ;
}

void StochasticGradient::display(){
std::cout << "Iterations : " << niterate << endl;
std::cout << "VaR-" << _alpha  << " : " << _xi << std::endl ;
std::cout << "CVaR-" << _alpha << " : " << _c << std::endl;
}

double StochasticGradient::H1(double fx){
	return (1.-(1./(1.-_alpha))*(fx > _xi));
}

double StochasticGradient::H2(double fx){
	return _c - _xi - (1./(1.-_alpha))*fmax(0., fx - _xi) ;
}







