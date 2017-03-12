#include "stochastic_process.hpp"

StochasticProcess::StochasticProcess(double T, double Nt, double X0)
{
	this->_T = T ;
	this->_Nt = Nt ;
	this->_x0 = X0 ;
	this->_h = T/Nt ;
}

