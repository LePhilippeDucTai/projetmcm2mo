#include "stochastic_gradient.hpp"

StochasticGradient::StochasticGradient(double gamma0,double (*gamma)(int,double), double xi0, double alpha) 
		: _gamma0(gamma0), _gamma(gamma), _xi(xi0), _c(0.), niterate(0.) {}

// double StochasticGradient::gamma(int n){
// 	return (this->_gamma0)/((double)n) ;
// }

// Cette fonction sera appelée dans une boucle pour i
// Elle calculera le taux de variation de la dernière valeur de xi avec la nouvelle calculée
// On s'arrêtera dans la boucle lorsque le taux d'erreur sera inférieur à une certaine valeur
void StochasticGradient::Iterate(int n, double fX){
	this->_lastxi = this->_xi ; // Enregistrement de la dernière valeur de xi
	this->_lastc = this->_c ; // Enregistrement de la dernière valeur de c
	double h1 = StochasticGradient::H1(fX);
	double h2 = StochasticGradient::H2(fX); 
	this->_xi = this->_xi - this->_gamma(n,this->_gamma0)*h1 ;
	this->_c = this->_c - this->_gamma(n,this->_gamma0)*h2 ;
	this->niterate++;
}

//Calcul de la différence |xi_(n+1)-xi_(n)| ;
double StochasticGradient::precision_xi(){
	return abs(this->_xi - this->_lastxi) ;
}

double StochasticGradient::precision_c(){
	return abs( this->_c - this->_lastc) ;
}

// Définition des accesseurs
double StochasticGradient::getVaR(){
	return this->_xi ;
}

double StochasticGradient::getCVaR(){
	return this->_c ;
}

void StochasticGradient::display(){
std::cout << "Iterations : " << this->niterate << endl;
std::cout << "VaR-" << this->_alpha  << " : " << this->_xi << std::endl ;
std::cout << "CVaR-" << this->_alpha << " : " << this->_c << std::endl;
}

double StochasticGradient::H1(double fx){
	return (1-(1/(1-this->_alpha))*(fx > this->_xi));
}

double StochasticGradient::H2(double fx){
	return this->_c - this->_xi - (1/(1-this->_alpha))*fmax(0, fx - this->_xi) ;
}








