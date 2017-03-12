#include "random_singleton.h"
#include <vector>
// Base class 
class StochasticProcess {
public :
    //Accesseurs
    virtual std::vector<double> get_X = 0 ;
    virtual std::vector<double> get_tgrid = 0 ;
    virtual get_h = 0 ;

private:
    double _T ; // Horizon de temps
    double _Nt ; // Niveau de discrétisation
    double _x0 ; //Valeur initiale
    double h ; // Pas de discrétisation
    std::vector <double> X ; // {X_{t_0}=X_0, ... , X_{t_n}} valeurs du processus discrétisé
    std::vector <double> tgrid ; // {0=t_0, t_1, ... , t_n = T} maillage du temps
}
