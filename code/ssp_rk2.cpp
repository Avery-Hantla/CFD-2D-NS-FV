/////////////////////////////////////////////////////////
//                  SSP RK2 Function
/////////////////////////////////////////////////////////
#include <vector>
#include <algorithm>

#include "class_q.hpp"
#include "class_mesh.hpp"
#include "class_flow.hpp"
#include "struct_size.hpp"
#include "struct_inputs.hpp"
#include "struct_BC.hpp"
#include "struct_residual.hpp"

#include "res.cpp"

void ssp_rk2(class_mesh* mesh, class_Q* Qbar, class_Q* Qface_c1, class_Q* Qface_c2, struct_residual* residual, struct_size* size, struct_inputs* inputs, class_flow* freestream, struct_BC* BC) {
    // Initilize Variables
    class_Q Q_star;
    Q_star.init(size->num_cells);

    // Compute Time Step
    double dt;
    dt = 0.001;

    // Calcualte Q star
    res(residual, mesh, size, inputs, Qbar, Qface_c1, Qface_c2, freestream, BC);
    for (int idx = 0; idx < size->num_cells; idx ++) {
        Q_star.p1[idx] = Qbar->p1[idx] + dt*residual->p1[idx];
        Q_star.p2[idx] = Qbar->p2[idx] + dt*residual->p2[idx];
        Q_star.p3[idx] = Qbar->p3[idx] + dt*residual->p3[idx];
        Q_star.p4[idx] = Qbar->p4[idx] + dt*residual->p4[idx];
    }

    // Calcualte Q n+1
    res(residual, mesh, size, inputs, &Q_star, Qface_c1, Qface_c2, freestream, BC);
    for (int idx = 0; idx < size->num_cells; idx ++) {
        Qbar->p1[idx] = Qbar->p1[idx] + Q_star.p1[idx] + dt*residual->p1[idx];
        Qbar->p2[idx] = Qbar->p2[idx] + Q_star.p2[idx] + dt*residual->p2[idx];
        Qbar->p3[idx] = Qbar->p3[idx] + Q_star.p3[idx] + dt*residual->p3[idx];
        Qbar->p4[idx] = Qbar->p4[idx] + Q_star.p4[idx] + dt*residual->p4[idx];
    }
}