#include "stochastic_gradient.hpp"
//Vectorized stochastic gradient

// StochasticGradient::StochasticGradient(double gamma0,double (*gamma)(int,double), double xi0, double alpha) 
// 		: _gamma0(gamma0), _gamma(gamma), _xi(xi0), _alpha(alpha), _c(0.), niterate(0.) {}

// Sequences are initialized at 0 by default
StochasticGradient::StochasticGradient(GaussianVector Y, double gamma0, double (*gamma)(int,double), double (*phi)(arma::vec), double alpha)
	: X(Y), _gamma0(gamma0), _gamma(gamma), _phi(phi),_alpha(alpha), _dim(Y.Size()),_c(0.),_xi(0.), 
		_theta(arma::zeros<arma::vec>(Y.Size())), _mu(arma::zeros<arma::vec>(Y.Size())), niterate(0){}

	// In case we need to initialize the sequences
	void StochasticGradient::init_seq(double xi0, double c0, arma::vec theta, arma::vec mu, double rho,double a,double b){
		_xi = xi0 ;
		_c = c0 ;
		_theta = theta ;
		_mu = mu ;
	}

// double StochasticGradient::gamma(int n){
// 	return (_gamma0)/((double)n) ;
// }

// Cette fonction sera appelée dans une boucle pour i
// Elle calculera le taux de variation de la dernière valeur de xi avec la nouvelle calculée
// On s'arrêtera dans la boucle lorsque le taux d'erreur sera inférieur à une certaine valeur
void StochasticGradient::Iterate(double epsilon){
	// Faire une boucle tant que la précision de xi + c est supérieure à epsilon
	int n = 1 ;
	double h1,h2 ;
	arma::vec simX ;
	do {
		_lastxi = _xi ; // Enregistrement de la dernière valeur de xi
		_lastc = _c ; // Enregistrement de la dernière valeur de c
		simX = X();
		h1 = StochasticGradient::H1(simX);
		h2 = StochasticGradient::H2(simX); 
		_xi = _xi - _gamma(n,_gamma0)*h1 ;
		_c = _c - _gamma(n,_gamma0)*h2 ;
		niterate++;
		n++ ;
	}while( precision_xi() + precision_c() > epsilon );
}


void StochasticGradient::IterateIS(double epsilon) {

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

double StochasticGradient::H1(arma::vec &x){
	return (1.-(1./(1.-_alpha))*(_phi(x) > _xi));
}

double StochasticGradient::H2(arma::vec &x){
	return _c - _xi - (1./(1.-_alpha))*fmax(0., _phi(x) - _xi) ;
}

// double Stochastic::Gradient::L1(arma::vec x){
//  return 0. ;
// }
// double Stochastic::Gradient::L2(arma::vec x){
//  return 0. ;
// }
// double Stochastic::Gradient::K1(arma::vec x){
//  return 0. ;
// }
// double Stochastic::Gradient::K2(arma::vec x){
//  return 0. ;
// }

// TODO LIST :
/*
- Inclure les constantes rho, a et b dans les RandomVect pour chaque loi.
- Ajouter le gradient de la pdf pour chaque loi dans RandomVect

- Dans le constructeur mettre en entrée un RandomVect : l'objet RandomVect contient toute l'information
nécessaire : génération des vecteurs aléatoires, contient la PDF, contient les constantes rho, a et b
nécessaires à l'algorithme d'Importance Sampling

- L'objet RandomVect faisant partie de l'algorithme, on peut réaliser la boucle de calculs dans la fonction Iterate

- Implémenter les fonctions L1, L2, K1, K2. Ces 4 fonctions prennent juste un argument de type arma::vec
Il n'est pas nécessaire d'ajouter d'autres champs car les valeurs de xi, C, mu et theta sont stockées dans la classe.


*/


