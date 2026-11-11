#include "euler.h"
#include "physic/physic.h"

cc::mat5 roe_abs_n(const cc::face_class& face){
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const double sl = std::sqrt(L.rho), sr = std::sqrt(R.rho);
    const double den = sl + sr;
    const double Hl = ph::Cp*L.T + 0.5*(L.u*L.u + L.v*L.v + L.w*L.w);
    const double Hr = ph::Cp*R.T + 0.5*(R.u*R.u + R.v*R.v + R.w*R.w);
    const double ubl = (sl*L.u + sr*R.u)/den;
    const double vbl = (sl*L.v + sr*R.v)/den;
    const double wbl = (sl*L.w + sr*R.w)/den;
    const double Hbl = (sl*Hl + sr*Hr)/den;
    const double abl = std::sqrt((ph::gamma - 1.0)*(Hbl - 0.5*(ubl*ubl + vbl*vbl + wbl*wbl)));
    cc::cell_class st;
    st.phy.u = ubl;
    st.phy.v = vbl;
    st.phy.w = wbl;
    st.phy.T = (Hbl - 0.5*(ubl*ubl + vbl*vbl + wbl*wbl))/ph::Cp;
    const cc::mat5 A = face.toface_jacobi(form_pfgpw_convect(st));
    const double len = face.nor.norm();
    const double un = ubl*face.nor.x + vbl*face.nor.y + wbl*face.nor.z;
    const double al = abl*len, ent = 0.1*abl*len;
    const double l1 = std::max(std::abs(un - al), ent);
    const double l4 = std::max(std::abs(un + al), ent);
    const double l0 = std::max(std::abs(un), ent);
    const double alpha = (l4 - l1)/(2.0*al);
    const double beta = (l4 + l1 - 2.0*l0)/(2.0*al*al);
    cc::mat5 M;
    for(int r = 0; r < 5; r++){
        for(int c = 0; c < 5; c++) M.a[r][c] = A.a[r][c] - (r == c ? un : 0.0);
    }
    cc::mat5 out;
    for(int r = 0; r < 5; r++){
        for(int c = 0; c < 5; c++){
            double m2 = 0.0;
            for(int k = 0; k < 5; k++) m2 += M.a[r][k]*M.a[k][c];
            out.a[r][c] = l0*(r == c ? 1.0 : 0.0) + alpha*M.a[r][c] + beta*m2;
        }
    }
    return out;
}