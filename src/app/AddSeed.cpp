#include "app/App.h"
#include "physics/Physics.h"

#include <cmath>

void add_seed(double amp){
    for(cc::cell_class& c : cc::CellList){
        const double r = std::hypot(c.center.x, c.center.y);
        if(r < 0.5 || r > 3.0) {continue;}
        const double wgt = amp*std::exp(-r*r/4.0);
        c.phy.u -= wgt*c.center.y;
        c.phy.v += wgt*c.center.x;
        *c.otphy = ph::form_otphy(c.phy);
        c.form_conservative();
    }
}
