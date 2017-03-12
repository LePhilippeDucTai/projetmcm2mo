#include "random_singleton.h"
#include <vector>
// Base class 
class StochasticProcess {
public :
    //Accesseurs
   	StochasticProcess(double T, double Nt, double X0) ;
    std::vector<double> get_X() ;
    std::vector<double> get_tgrid()  ;
    double get_h() ;

private:
    double _T ; // Horizon de temps
    double _Nt ; // Niveau de discrétisation
    double _x0 ; //Valeur initiale
    double _h ; // Pas de discrétisation
    std::vector <double> X ; // {X_{t_0}=X_0, ... , X_{t_n}} valeurs du processus discrétisé
    std::vector <double> tgrid ; // {0=t_0, t_1, ... , t_n = T} maillage du temps

};
