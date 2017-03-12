// random-singleton.cpp (Singleton)
// Implementation de la classe Random / Implementation of the class Random
#include "random_singleton.h"

#include <cfloat> //DBL_EPSILON
using namespace std;

inline static double sqr(double x) {return x*x;}

//Ces constantes sont requises par l'algorithme du generateur aleatoire
//These consts are required by the algorithm
const int Random::NTAB = 32;

const long Random::IM1 = 2147483563;
const long Random::IM2 = 2147483399;
const long Random::IMM1 = Random::IM1 - 1;
const double Random::AM = 1.0/Random::IM1;
const int Random::IA1 = 40014;
const int Random::IA2 = 40692;
const int Random::IQ1 = 53668;
const int Random::IQ2 = 52774;
const int Random::IR1 = 12211;
const int Random::IR2 = 3791;
const int Random::NDIV = 1 + Random::IMM1/Random::NTAB;
const double Random::EPS = DBL_EPSILON;
const double Random::RNMX = 1.0-Random::EPS;

long  Random::idum   = 0;
long  Random::idum2  = 123456789L;
long  Random::iy     = 0;
long* Random::iv     = 0;

//definition du singleton
Random Random::Singleton(0);

void Random::Randomize(long seed)
{
  idum = (seed <= 0) ? (seed == 0 ? 1 : -seed) : seed; //be sure to prevent idum=0

  idum2 = idum;
  if (idum2 < 0) idum2 = -idum2;
   
  long k = 0;
  for(int j = NTAB+7 ; j >= 0 ; j--)
  {
    k = idum/IQ1;
    idum = IA1*(idum - k*IQ1) - k*IR1;
    if (idum < 0)
      idum += IM1;

    if (j < NTAB)
      iv[j] = idum;
  }

  iy = iv[0];
}
//fin Randomize()

//Cette fonction renvoie un double pseudo-aleatoire uniformement dans ]0;1[
//Il s'agit de l'algorithme de L'Ecuyer avec melange de Bays-Durham
//This function returns a double, taken uniformly in ]0;1[
//It is the algorithm of L'Ecuyer with a Bays-Durham shuffle
double Random::theRandom(void)
{
  long k = idum/IQ1;
  idum = IA1*(idum - k*IQ1) - k*IR1;
  if (idum < 0)
    idum +=IM1;

  k = idum2/IQ2;
  idum2 = IA2*(idum2 -k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 +=IM2;

  long j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy < 1)
    iy += IMM1;
      
  double temp = AM*iy;
  if (temp >= RNMX)
    return RNMX;  //empeche de renvoyer 1 / prevents from returning 1
  else
    return temp;
}
//fin theRandom()

//Renvoie un double de la Gaussienne de moyenne et d'ecartype specifies
//Returns a double taken on a Gaussian with specified mean and standard deviation
double Random::Gaussian(double mean, double standardDeviation){
  double x = Random::Uniform();
  int signe = 0;
  double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328, d1 = 1.432788, d2 = 0.189269, d3 = 0.001308 ;
  if(x>0.5){
    signe = 1 ;
    x=1.0-x;
  }
  else 
    signe=-1.0;
  double t = sqrt(-2.0*log(x));
  return mean + standardDeviation*(signe * (t-((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t)))));
}

static double inverse_N(double x) 
{
int signe = 0 ;
  double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328, d1 = 1.432788, d2 = 0.189269, d3 = 0.001308 ;
  if(x>0.5){
    signe = 1 ;
    x=1.0-x;
  }
  else 
    signe=-1.0;
  double t = sqrt(-2.0*log(x));
  return(signe * (t-((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t)))));
}

double Random::LogNormal(double mu, double sigma)
{
  return( exp(Random::Gaussian(mu,sigma))) ;
}

double Random::Exponential(double lambda)
{
  return (-1)*log(Random::Uniform())/lambda ;
}




//fin random-singleton.cpp
