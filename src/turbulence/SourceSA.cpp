#include "physics/Physics.h"
#include "turbulence/Turbulence.h"

double SA::source_SA(cc::cell_class& cell){
    build_SA(cell);
    const double d2 = cell.tur->sad*cell.tur->sad;
    const double production = Cb1*(1.0 - cell.tur->ft2)*cell.tur->Sbl*cell.phy.rho*cell.tur->miubl;
    const double destruction = cell.phy.rho*(Cw1*cell.tur->fw - Cb1/(kappa*kappa)*cell.tur->ft2)*
                               cell.tur->miubl*cell.tur->miubl*(1.0/d2);
    const double compressible = C5*cell.phy.rho*cell.tur->miubl*cell.tur->miubl*cell.tur->S2/
                                (ph::gamma*ph::R*cell.phy.T);
    return production - destruction - compressible;
}
