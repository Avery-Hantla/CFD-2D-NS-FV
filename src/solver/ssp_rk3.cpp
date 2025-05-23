/////////////////////////////////////////////////////////
//                  SSP RK2 Function
/////////////////////////////////////////////////////////
#include <vector>
#include <algorithm>

#include "../support/class_q.hpp"
#include "../support/class_mesh.hpp"
#include "../support/class_flow.hpp"
#include "../support/class_residual.hpp"
#include "../support/struct_size.hpp"
#include "../support/struct_inputs.hpp"
#include "../support/struct_BC.hpp"
#include "../support/struct_time.hpp"

#include "../include/res.hpp"

void ssp_rk3(class_mesh* mesh, class_Q* Qbar, class_Q* Qface, class_Q* Qface_c1, class_Q* Qface_c2, class_residual* residual, struct_size* size, struct_inputs* inputs, class_flow* freestream, struct_BC* BC, struct_time* time) {
    // Initilize Variables
    class_Q Q_star;
    Q_star.init(size->num_cells, freestream->gamma, freestream->R);

    // Compute Time Step
    int c1, c2;
    std::vector<double> dt;
    dt.assign(size->num_cells, 0);
    double temp;
    for (int idx = 0; idx < size->num_faces; idx ++) {
        c1 = mesh->face_cell1[idx];
        c2 = mesh->face_cell2[idx];

        dt[c1] += (std::abs(Qface_c1->Vn_avg[idx])+Qface_c1->c_avg[idx])*mesh->face_area[idx];
        if (c2 >= 0){
            dt[c2] += (std::abs(Qface_c2->Vn_avg[idx])+Qface_c2->c_avg[idx])*mesh->face_area[idx];
        }
    }
    
    if (time->type == 1) { // GLOBAL Time step
        double dt_temp = *std::max_element(dt.begin(),dt.end());
        for (int idx = 0; idx < size->num_cells; idx ++) {
            dt[idx] = (2*time->CFL*mesh->cell_vol[idx])/dt_temp;
        }
    } else if (time->type == 2) { // LOCAL Time step
        for (int idx = 0; idx < size->num_cells; idx ++) {
            dt[idx] = (2*time->CFL*mesh->cell_vol[idx])/dt[idx];
        }
    }

    // Calcualte Q star 1
    res(residual, mesh, size, inputs, Qbar, Qface, Qface_c1, Qface_c2, freestream, BC);
    for (int idx = 0; idx < size->num_cells; idx ++) {
        Q_star.p1[idx] = Qbar->p1[idx] + dt[idx]*residual->p1[idx];
        Q_star.p2[idx] = Qbar->p2[idx] + dt[idx]*residual->p2[idx];
        Q_star.p3[idx] = Qbar->p3[idx] + dt[idx]*residual->p3[idx];
        Q_star.p4[idx] = Qbar->p4[idx] + dt[idx]*residual->p4[idx];
    }
    Q_star.updateflow();

    // Calcualte Q star 2
    res(residual, mesh, size, inputs, &Q_star, Qface, Qface_c1, Qface_c2, freestream, BC);
    for (int idx = 0; idx < size->num_cells; idx ++) {
        Q_star.p1[idx] = 0.75*Qbar->p1[idx] + 0.25*Q_star.p1[idx] + 0.25*dt[idx]*residual->p1[idx];
        Q_star.p2[idx] = 0.75*Qbar->p2[idx] + 0.25*Q_star.p2[idx] + 0.25*dt[idx]*residual->p2[idx];
        Q_star.p3[idx] = 0.75*Qbar->p3[idx] + 0.25*Q_star.p3[idx] + 0.25*dt[idx]*residual->p3[idx];
        Q_star.p4[idx] = 0.75*Qbar->p4[idx] + 0.25*Q_star.p4[idx] + 0.25*dt[idx]*residual->p4[idx];
    }
    Q_star.updateflow();

    // Calcualte Q n+1
    res(residual, mesh, size, inputs, &Q_star, Qface, Qface_c1, Qface_c2, freestream, BC);
    for (int idx = 0; idx < size->num_cells; idx ++) {
        Qbar->p1[idx] = (1.0/3.0)*Qbar->p1[idx] + (2.0/3.0)*Q_star.p1[idx] + (2.0/3.0)*dt[idx]*residual->p1[idx];
        Qbar->p2[idx] = (1.0/3.0)*Qbar->p2[idx] + (2.0/3.0)*Q_star.p2[idx] + (2.0/3.0)*dt[idx]*residual->p2[idx];
        Qbar->p3[idx] = (1.0/3.0)*Qbar->p3[idx] + (2.0/3.0)*Q_star.p3[idx] + (2.0/3.0)*dt[idx]*residual->p3[idx];
        Qbar->p4[idx] = (1.0/3.0)*Qbar->p4[idx] + (2.0/3.0)*Q_star.p4[idx] + (2.0/3.0)*dt[idx]*residual->p4[idx];
    }
    Qbar->updateflow();
}