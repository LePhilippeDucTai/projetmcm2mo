#include "stochastic_gradient.hpp"
//Vectorized stochastic gradient

// StochasticGradient::StochasticGradient(double gamma0,double (*gamma)(int,double), double xi0, double alpha) 
// 		: _gamma0(gamma0), _gamma(gamma), _xi(xi0), _alpha(alpha), _c(0.), niterate(0.) {}

// Sequences are initialized at 0 by default
StochasticGradient::StochasticGradient( GaussianVector Y, 
																				double gamma0, 
																				double (*gamma)(int,double), 
																				double (*phi)(const arma::vec &),
																				double alpha)
				: X(Y), _gamma0(gamma0), 
					_gamma(gamma), _phi(phi), 
					_alpha(alpha), _dim(Y.Size()),
					_c(0.),_xi(0.), _theta(arma::zeros<arma::vec>(Y.Size())), 
					_mu(arma::zeros<arma::vec>(Y.Size())), 
					niterate(0) ,
					_cesaro_xi(0.),
					_cesaro_c(0.),
					_varemp_xi(0.),
					_varemp_c(0.)
					{}

	// In case we need to initialize the sequences
	void StochasticGradient::init_seq(double xi0, double c0, 
																		arma::vec theta, arma::vec mu, 
																		double rho,double a,double b) {
		_xi = xi0 ;
		_c = c0 ;
		_theta = theta ;
		_mu = mu ;
	}



void StochasticGradient::Calculation(double &h1, double &h2, arma::vec &simX){
	_lastxi = _xi ; // Enregistrement de la dernière valeur de xi
	_lastc = _c ; // Enregistrement de la dernière valeur de c
	simX = X();
	niterate++;
	h1 = StochasticGradient::H1(simX);
	h2 = StochasticGradient::H2(simX); 
	_xi = _xi - _gamma(niterate,_gamma0)*h1 ;
	_c = _c - _gamma(niterate,_gamma0)*h2 ;
	_cesaro_xi = _cesaro_xi + _xi ;
	_cesaro_c = _cesaro_c + _c ;
	_varemp_xi = _varemp_xi + (_xi - _cesaro_xi)*(_xi - _cesaro_xi) ;
	_varemp_c = _varemp_c + (_c - _cesaro_c)*(_c - _cesaro_c);
}

void StochasticGradient::CalculationIS(double &l1, double &l2, arma::vec &k1, arma::vec &k2, arma::vec &simX){
	_lastxi = _xi ; // Enregistrement de la dernière valeur de xi
	_lastc = _c ; // Enregistrement de la dernière valeur de c
	simX = X();
	niterate++;
	l1 = L1(simX);
	l2 = L2(simX);
	k1 = K1(simX);
	k2 = K2(simX);
	_xi = _xi - _gamma(niterate,_gamma0)*l1 ;
	_c = _c - _gamma(niterate,_gamma0)*l2 ;
	_cesaro_xi = _cesaro_xi + _xi ;
	_cesaro_c = _cesaro_c + _c ;
	_theta = _theta - _gamma(niterate,_gamma0)*k1 ;
	_mu = _mu - _gamma(niterate,_gamma0)*k2 ;
	_varemp_xi = _varemp_xi + (_xi - _cesaro_xi)*(_xi - _cesaro_xi) ;
	_varemp_c = _varemp_c + (_c - _cesaro_c)*(_c - _cesaro_c);
}

void StochasticGradient::EmpiricalMeans()
{
	_cesaro_xi = _cesaro_xi/niterate ;
	_cesaro_c = _cesaro_c/niterate ;
	_varemp_xi = _varemp_xi/niterate ;
	_varemp_c = _varemp_c/niterate ;
}


// Cette fonction sera appelée dans une boucle pour i
// Elle calculera le taux de variation de la dernière valeur de xi avec la nouvelle calculée
// On s'arrêtera dans la boucle lorsque le taux d'erreur sera inférieur à une certaine valeur
void StochasticGradient::Iterate(double epsilon){
	// Faire une boucle tant que la précision de xi + c est supérieure à epsilon
	int n = 1 ;
	double h1,h2 ;
	arma::vec simX ;
	do {
		Calculation(h1,h2,simX);
	}while( precision_xi() + precision_c() > epsilon );
	EmpiricalMeans() ;
}

void StochasticGradient::Iterate(int nsim){
	// Faire une boucle tant que la précision de xi + c est supérieure à epsilon
	int n = 1 ;
	double h1,h2 ;
	arma::vec simX ;
	do {
		Calculation(h1,h2,simX);
	}while( niterate < nsim );
	EmpiricalMeans() ;
}

void StochasticGradient::IterateIS(double epsilon) {
	// Faire une boucle tant que la précision de xi + c est supérieure à epsilon
	double l1,l2 ;
	arma::vec k1,k2;
	arma::vec simX ;
	do {
		CalculationIS(l1,l2,k1,k2,simX);
	}while( precision_xi() + precision_c() > epsilon );
	EmpiricalMeans() ;
}
void StochasticGradient::IterateIS(int nsim) {
	// Faire une boucle tant que la précision de xi + c est supérieure à epsilon
	double l1,l2 ;
	arma::vec k1,k2;
	arma::vec simX ;
	do {
		CalculationIS(l1,l2,k1,k2,simX);
	}while( niterate < nsim );
	EmpiricalMeans() ;
}


//Calcul de la différence |xi_(n+1)-xi_(n)| ;
inline double StochasticGradient::precision_xi() const{
	return abs(_xi - _lastxi) ;
}

inline double StochasticGradient::precision_c() const {
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
std::cout << "Cesaro VaR-" << _alpha << " : " << _cesaro_xi << std::endl ;
std::cout << "Cesaro CVaR-" << _alpha << " : " << _cesaro_c << std::endl ;
std::cout << "Variance Empirique VaR-" << _alpha << " : " << _varemp_xi << std::endl ;
std::cout << "Variance Empirique CVaR-" << _alpha << " : " << _varemp_c << std::endl;

}

double StochasticGradient::H1(const arma::vec &x){
	return (1.-(1./(1.-_alpha))*(_phi(x) > _xi));
}

double StochasticGradient::H2(const arma::vec &x){
	return _c - _xi - (1./(1.-_alpha))*fmax(0., _phi(x) - _xi) ;
}

double StochasticGradient::L1(const arma::vec &x){
double cste = exp(-1.*X.SGIS_rho()*pow(arma::norm(_theta),X.SGIS_b()));
arma::vec xpt = x+_theta ;
return cste*(1. - (1./(1.-_alpha))*(_phi(xpt)>= _xi)*X.Pdf(xpt)/X.Pdf(x));
}


double StochasticGradient::L2(const arma::vec &x){
 arma::vec xpm = x+_mu ;
 return _c - (_xi + (1./(1.-_alpha))*fmax(0.,_phi(xpm)-_xi)*X.Pdf(xpm)/X.Pdf(x));
}

// K1 is needed for theta, hence is a vector
arma::vec StochasticGradient::K1(arma::vec &x){
 double p = X.Pdf(x); // p(x)
 double pmt = X.Pdf(x-_theta); // p(x-theta)
 double pm2t = X.Pdf(x-2.*_theta); //p(x-2theta)
 arma::vec gradpm2t = X.Gradient(x-2.*_theta); //Grad.p(x-2theta)
 double cste = exp(-2.*X.SGIS_rho()*pow(arma::norm(_theta),X.SGIS_b()));
 arma::vec xmt = x-_theta ;
 arma::vec W1t = (_phi(xmt)>= _xi)*(cste*pmt*pmt/(p*pm2t*pm2t))*gradpm2t;  // W1 tilde
 return W1t ;
}

// K2 is needed for mu, hence is a vector
arma::vec StochasticGradient::K2(arma::vec &x){
 double p = X.Pdf(x); // p(x)
 double pmm = X.Pdf(x-_mu); // p(x-mu)
 double pm2m = X.Pdf(x-2.*_mu); //p(x-2mu)
 arma::vec gradpm2m = X.Gradient(x-2.*_mu); //Grad.p(x-2mu)
 double cste = exp(-2.*X.SGIS_rho()*pow(arma::norm(_theta),X.SGIS_b())-X.SGIS_a()*(pow(arma::norm(_mu),2)+1.));
 arma::vec xmm = x-_mu ; 
 arma::vec W2t = pow(fmax(0.,_phi(xmm)- _xi),2)*(cste*pmm*pmm/(p*pm2m*pm2m))*gradpm2m;  // W2tilde
 return W2t ;
}









