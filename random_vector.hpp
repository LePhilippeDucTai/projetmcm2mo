#ifndef RANDOM_VECTOR_HPP
#define RANDOM_VECTOR_HPP
#include <functional>
#include <vector>
#include <cmath>
#include <list>
#include <iostream>
#include <random>
#include <chrono>
#include "random_variable.hpp"

template <class T> struct RandomVect {
	
	typedef T result_type;
	RandomVect() 
			: value(0) {}
	RandomVect(T value)
			: value(value) {}
	virtual ~RandomVect() {}
	virtual T operator()() = 0; 
	virtual T Pdf(double) = 0 ;
	T Current() const {
		return value;
	}
 protected:
	T value;
	std::mt19937 generator;
};


struct UniformVector : public RandomVect<std::vector<double>> {
	UniformVector(double left_boundary , double right_boundary)
			: left_boundary(left_boundary), right_boundary(right_boundary), generator(seed) {}

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
	std::mt19937 generator;
};
















