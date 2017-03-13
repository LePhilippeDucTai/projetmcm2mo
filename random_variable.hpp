#ifndef RANDOM_VARIABLE_HPP
#define RANDOM_VARIABLE_HPP
#include <functional>
#include <vector>
#include <cmath>
#include <list>
#include <iostream>
#include <random>
#include <chrono>

static unsigned seed;
//TODO :
//	- Change all argument into vectors
//	- Do the pdf of multivariate normal distribution and also its derivative
template <class T> struct RandomVar {
	
	typedef T result_type;
	RandomVar() 
			: value(0) {}
	RandomVar(T value)
			: value(value),generator(1) {} // default seed = 1
	virtual ~RandomVar() {}
	virtual T operator()() = 0; 
	virtual T Pdf(double) = 0 ;
	T Current() const {
		return value;
	}

	void init(int seed_) {
		seed = seed_ ;
		generator.seed(seed);
	}

 protected:
	T value;
	int seed ;
	std::mt19937 generator;

};

struct Uniform : public RandomVar<double> {
	Uniform(double left_boundary , double right_boundary)
			: left_boundary(left_boundary), right_boundary(right_boundary) {}

	double operator()() {
		return value = left_boundary + (right_boundary - left_boundary)*
					(generator()/static_cast<double> (generator.max()));
	}
	inline double Pdf(double x) {
		return (x >= left_boundary && x <= right_boundary) ? 
							1./static_cast<double>(right_boundary - left_boundary) : 0.;
	}
 private:
	double left_boundary, right_boundary;

};

struct Gaussian : public RandomVar<double> {
	Gaussian(double mean, double std) 
			: mean(mean), std(std) {}
	double operator()() {
		double x = generator()/static_cast<double> (generator.max());
		int signe = 0;
		double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328, d1 = 1.432788, d2 = 0.189269, d3 = 0.001308 ;

		if(x>0.5){
			signe = 1 ;
			x=1.0-x;
		}
		else 
			signe=-1.0;
			
		double t = sqrt(-2.0*log(x));
		return mean + std*(signe * (t-((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t)))));
	}
	inline double Pdf(double x) {
		return exp(-0.5*(x-mean)*(x-mean)/std)/sqrt(2*M_PI*std);
	}
 private :
 double mean, std;

};

struct Exponential : public RandomVar<double> {
	Exponential(double lambda) 
			: lambda(lambda), inv_lambda(1./lambda){}
	double operator()() {
		return -inv_lambda*log(generator()/static_cast<double> (generator.max()));
	}
		inline double Pdf(double x) {
		return (x > 0) ? lambda*exp(-lambda*x): 0;
	}
 private :
	double inv_lambda, lambda;

};

struct LogNormal : public RandomVar<double> {
	LogNormal(double mean, double std)
		: mean(mean), std(std), G(mean,std) {}
	double operator()() {
		return exp(G());
	}
		inline double Pdf(double x) {
		return (x > 0) ? exp(-0.5*(log(x)-mean)*(log(x)-mean)/std)/(sqrt(2*M_PI*std)*x): 0;
	}
 private :
 	Gaussian G;
 	double mean;
 	double std;
};

struct ChiSquared : public RandomVar<double> {
	ChiSquared(int N) : N(N), G(0,1) {};
	double operator()() {
		value = 0;
		for (int j = 0; j < N; j++) 
			value += G()*G.Current(); // This allows us to take the Current gaussian
		return value;
	};
//TODO : Do the pdf of Chi-Squared distribution
		inline double Pdf(double x) {
		return (x > 0) ? 1: 0;
	}
	private:
		int N;
		Gaussian G;
};

struct Gamma : public RandomVar<double> {
	Gamma(double a, double b) 
		: a(a), b(b), Z(0,1), U(0,1) {}
	double operator()() {

  if( a < 1) {
  	Gamma G(a+1,b);
 		return G()*pow(U(),1./a);
 	}
 	double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt (d);
  double z = Z();
  double u = U();
  double v = pow((1+c*z),3.);
	while (	z < -1./c && (log(u) > 0.5*z*z + d*v + d*log(v))) {
		z = Z();
		u = U();
		v = pow((1+c*z),3.);
	}

	return d*v/b;
}	
//TODO : Do the pdf of Gamma distribution
		inline double Pdf(double x) {
		return (x > 0) ? 1: 0;
	}
 private :
 double a, b;
 Gaussian Z;
 Uniform U;
};
		


template <typename Gen>
std::vector<double> MonteCarlo(int n, Gen G, double (*phi)(double))
{
	std::vector<double> result(3,0);
	double x;
	for (int j = 0; j < n; j++) {
		x = G();
		phix = phi(x) ;
		result[0] += phix;
		result[1] += phix*phix;
	}
	result[0] /= (double) n;
	result[1] = (result[1] - n*result[0]*result[0])/(double)(n-1);
	result[2] = 1.96*sqrt(result[1]/(double) n);
	return result;
}


//	void gen_Poisson(int nsim, double lambda) ;
//	void gen_Laplace(int nsim,double lambda) ;

/*ReturnType LongClassName::ReallyReallyReallyLongFunctionName(
		Type par_name1,  // 4 space indent
    Type par_name2,
    Type par_name3) {
  DoSomething();  // 2 space indent
  ...
}

class MyClass : public OtherClass {
 public:      // Note the 1 space indent!
  MyClass();  // Regular 2 space indent.
  explicit MyClass(int var);
  ~MyClass() {}

  void SomeFunction();
  void SomeFunctionThatDoesNothing() {
  }

  void set_some_var(int var) { some_var_ = var; }
  int some_var() const { return some_var_; }

 private:
  bool SomeInternalFunction();

  int some_var_;
  int some_other_var_;
};

*/














#endif
