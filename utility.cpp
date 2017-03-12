#include "utility.hpp"


std::vector<double> range(double a, double b, double by){
	std::vector<double> res ;
	double start ;
	start = a ;
	while(start < b) 
	{
		res.push_back(start);
		start = a + (b-a)/(std::abs(b-a))*by ;
	}
	return res ;
}

// template<typename T> 
// void display_vect(std::vector<T> t){
// 	for( std::vector<double>::iterator it =t.begin() ; it != t.end() ; ++it ){
// 		std::cout << *it << ' ' ;
// 	}
// 	std::cout << std::endl;
// }

template<typename T> 
void display_vect(std::vector<T> t){
    for(auto i : t){
       std::cout << i ;
    }
    std::cout << "\n";
}



