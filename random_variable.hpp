#include "random_singleton.h"

/*
La classe RandomVariable permet d'assigner à un objet X de type RandomVariable
une certaine loi spécifique. 
*/
class RandomVariable
{
public:
	RandomVariable::RandomVariable(int seed);
	void gen_Normal(int nsim, double mu, double sigma) ;
	void gen_Uniform(int nsim, double a, double b) ;
	void gen_Exponential(int nsim, double lambda) ;
	void gen_LogNormal(int nsim, double mu, double sigma) ;
	void gen_Poisson(int nsim, double lambda);
	void gen_Laplace(int nsim,double lambda);
	void MonteCarlo( double(*f)(std::vector<double>));
	

private:
	std::vector<double> value();
	double E ; //Esperance empirique
	double V ; //Variance empirique

}
