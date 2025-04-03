/////////////////////////////////////////////////////////
//           Function To Compute Roe Flux
/////////////////////////////////////////////////////////

#include "class_q.hpp"
#include "class_mesh.hpp"
#include "class_f.hpp"
#include "class_flow.hpp"
#include "struct_size.hpp"
#include "struct_BC.hpp"

void roe(class_F* F_roe, class_Q* Qface_c1, class_Q* Qface_c2, class_mesh* mesh, class_flow* freestream, struct_size* size, struct_BC* BC) {
    class_F F_c1, F_c2;
    F_c1.init(size->num_faces, freestream->gamma);
    F_c2.init(size->num_faces, freestream->gamma);
    
    F_c1.update(mesh, Qface_c1, size); 
    F_c2.update(mesh, Qface_c2, size);

    for (int idx = 0; idx < size->num_faces; idx++) {
        // Find Vn and c averag 
        double Vn_avg = (F_c1.Vn[idx] + F_c2.Vn[idx])/2;
        double c_avg = (Qface_c1->c[idx] + Qface_c2->c[idx])/2;

        Qface_c1->Vn_avg[idx] = Vn_avg;
        Qface_c1->c_avg[idx] = c_avg;
        Qface_c2->Vn_avg[idx] = Vn_avg;
        Qface_c2->c_avg[idx] = c_avg;

        // Compute Rusanov Flux
        F_roe->p1[idx] = (F_c1.p1[idx]+F_c2.p1[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p1[idx] - Qface_c1->p1[idx]);
        F_roe->p2[idx] = (F_c1.p2[idx]+F_c2.p2[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p2[idx] - Qface_c1->p2[idx]);
        F_roe->p3[idx] = (F_c1.p3[idx]+F_c2.p3[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p3[idx] - Qface_c1->p3[idx]);
        F_roe->p4[idx] = (F_c1.p4[idx]+F_c2.p4[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p4[idx] - Qface_c1->p4[idx]);
    }
}