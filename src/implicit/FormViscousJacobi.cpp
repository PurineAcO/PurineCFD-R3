#include "implicit/Implicit.h"
#include "euler/Euler.h"
#include "turbulence/Turbulence.h"

void impl::form_viscous_jacobi(cc::face_class& face){
    face.imats.Bf = SA::diffusion_state_jac(face);
    if(face.type != cc::INTER){
        face.imats.Ab.clear();
        face.imats.Tf.clear();
        return;
    }
    face.imats.Ab = roe_abs_n(face);
    face.imats.Tf = SA::diffusion_grad_jac(face, SA::normal_dist(&face));
}
