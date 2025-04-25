#include <vector>
#include <algorithm>
// #include <cmath>

#include "class_mesh.hpp"
#include "class_q.hpp"
#include "struct_size.hpp"

void squeeze(class_mesh* mesh, class_Q* Qbar, class_Q* Qface_c1, class_Q* Qface_c2, struct_size* size) {
    std::vector<double> Qmin_p1(size->num_cells, 100000000000), Qmin_p2(size->num_cells, 100000000000), Qmin_p3(size->num_cells, 100000000000), Qmin_p4(size->num_cells, 100000000000);
    std::vector<double> Qmax_p1(size->num_cells, -100000000000), Qmax_p2(size->num_cells, -100000000000), Qmax_p3(size->num_cells, -100000000000), Qmax_p4(size->num_cells, -100000000000);
    
    double temp_Qmin_p1, temp_Qmin_p2, temp_Qmin_p3, temp_Qmin_p4;
    double temp_Qmax_p1, temp_Qmax_p2, temp_Qmax_p3, temp_Qmax_p4;

    // Find Qmax and Qmin for cell
    for (int idx = 0; idx < size->num_faces; idx ++) {
        int cell1 = mesh->face_cell1[idx];
        int cell2 = mesh->face_cell2[idx];

        if (cell2 >= 0) {
            temp_Qmin_p1 = std::min({ Qbar->p1[cell1], Qbar->p1[cell2] });
            temp_Qmin_p2 = std::min({ Qbar->p2[cell1], Qbar->p2[cell2] });
            temp_Qmin_p3 = std::min({ Qbar->p3[cell1], Qbar->p3[cell2] });
            temp_Qmin_p4 = std::min({ Qbar->p4[cell1], Qbar->p4[cell2] });

            temp_Qmax_p1 = std::max({ Qbar->p1[cell1], Qbar->p1[cell2] });
            temp_Qmax_p2 = std::max({ Qbar->p2[cell1], Qbar->p2[cell2] });
            temp_Qmax_p3 = std::max({ Qbar->p3[cell1], Qbar->p3[cell2] });
            temp_Qmax_p4 = std::max({ Qbar->p4[cell1], Qbar->p4[cell2] });

            Qmin_p1[cell2] = std::min({ Qmin_p1[cell2], temp_Qmin_p1 });
            Qmin_p2[cell2] = std::min({ Qmin_p2[cell2], temp_Qmin_p2 });
            Qmin_p3[cell2] = std::min({ Qmin_p3[cell2], temp_Qmin_p3 });
            Qmin_p4[cell2] = std::min({ Qmin_p4[cell2], temp_Qmin_p4 });

            Qmax_p1[cell2] = std::max({ Qmax_p1[cell2], temp_Qmax_p1 });
            Qmax_p2[cell2] = std::max({ Qmax_p2[cell2], temp_Qmax_p2 });
            Qmax_p3[cell2] = std::max({ Qmax_p3[cell2], temp_Qmax_p3 });
            Qmax_p4[cell2] = std::max({ Qmax_p4[cell2], temp_Qmax_p4 });
        } else {
            temp_Qmin_p1 = Qbar->p1[cell1];
            temp_Qmin_p2 = Qbar->p2[cell1];
            temp_Qmin_p3 = Qbar->p3[cell1];
            temp_Qmin_p4 = Qbar->p4[cell1];
            
            temp_Qmax_p1 = Qbar->p1[cell1];
            temp_Qmax_p2 = Qbar->p2[cell1];
            temp_Qmax_p3 = Qbar->p3[cell1];
            temp_Qmax_p4 = Qbar->p4[cell1];
        }

        Qmin_p1[cell1] = std::min({ Qmin_p1[cell1], temp_Qmin_p1 });
        Qmin_p2[cell1] = std::min({ Qmin_p2[cell1], temp_Qmin_p2 });
        Qmin_p3[cell1] = std::min({ Qmin_p3[cell1], temp_Qmin_p3 });
        Qmin_p4[cell1] = std::min({ Qmin_p4[cell1], temp_Qmin_p4 });

        Qmax_p1[cell1] = std::max({ Qmax_p1[cell1], temp_Qmax_p1 });
        Qmax_p2[cell1] = std::max({ Qmax_p2[cell1], temp_Qmax_p2 });
        Qmax_p3[cell1] = std::max({ Qmax_p3[cell1], temp_Qmax_p3 });
        Qmax_p4[cell1] = std::max({ Qmax_p4[cell1], temp_Qmax_p4 });
    }

    // Find phi 
    std::vector<double> phi_p1_c1(size->num_faces,1), phi_p2_c1(size->num_faces,1), phi_p3_c1(size->num_faces,1),phi_p4_c1(size->num_faces,1);
    std::vector<double> phi_p1_c2(size->num_faces,1), phi_p2_c2(size->num_faces,1), phi_p3_c2(size->num_faces,1),phi_p4_c2(size->num_faces,1);

    for (int idx = 0; idx < size->num_faces; idx ++) {
        int cell1 = mesh->face_cell1[idx];
        int cell2 = mesh->face_cell2[idx];

        // Cell 1
        if (Qface_c1->p1[idx] < Qmin_p1[cell1]) {
            phi_p1_c1[idx] = (Qmin_p1[cell1] - Qbar->p1[cell1]) / (Qface_c1->p1[idx] - Qbar->p1[cell1]);
        } else if (Qface_c1->p1[idx] > Qmax_p1[cell1]) {
            phi_p1_c1[idx] = (Qmax_p1[cell1] - Qbar->p1[cell1]) / (Qface_c1->p1[idx] - Qbar->p1[cell1]);
        } 
        
        if (Qface_c1->p2[idx] < Qmin_p2[cell1]) {
            phi_p2_c1[idx] = (Qmin_p2[cell1] - Qbar->p2[cell1]) / (Qface_c1->p2[idx] - Qbar->p2[cell1]);
        } else if (Qface_c1->p2[idx] > Qmax_p2[cell1]) {
            phi_p2_c1[idx] = (Qmax_p2[cell1] - Qbar->p2[cell1]) / (Qface_c1->p2[idx] - Qbar->p2[cell1]);
        }

        if (Qface_c1->p3[idx] < Qmin_p3[cell1]) {
            phi_p3_c1[idx] = (Qmin_p3[cell1] - Qbar->p3[cell1]) / (Qface_c1->p3[idx] - Qbar->p3[cell1]);
        } else if (Qface_c1->p3[idx] > Qmax_p3[cell1]) {
            phi_p3_c1[idx] = (Qmax_p3[cell1] - Qbar->p3[cell1]) / (Qface_c1->p3[idx] - Qbar->p3[cell1]);
        }

        if (Qface_c1->p4[idx] < Qmin_p4[cell1]) {
            phi_p4_c1[idx] = (Qmin_p4[cell1] - Qbar->p4[cell1]) / (Qface_c1->p4[idx] - Qbar->p4[cell1]);
        } else if (Qface_c1->p4[idx] > Qmax_p4[cell1]) {
            phi_p4_c1[idx] = (Qmax_p4[cell1] - Qbar->p4[cell1]) / (Qface_c1->p4[idx] - Qbar->p4[cell1]);
        }

        // Cell 2
        if (Qface_c2->p1[idx] < Qmin_p1[cell2]) {
            phi_p1_c2[idx] = (Qmin_p1[cell2] - Qbar->p1[cell2]) / (Qface_c2->p1[idx] - Qbar->p1[cell2]);
        } else if (Qface_c2->p1[idx] > Qmax_p1[cell2]) {
            phi_p1_c2[idx] = (Qmax_p1[cell2] - Qbar->p1[cell2]) / (Qface_c2->p1[idx] - Qbar->p1[cell2]);
        } 

        if (Qface_c2->p2[idx] < Qmin_p2[cell2]) {
            phi_p2_c2[idx] = (Qmin_p2[cell2] - Qbar->p2[cell2]) / (Qface_c2->p2[idx] - Qbar->p2[cell2]);
        } else if (Qface_c2->p2[idx] > Qmax_p2[cell2]) {
            phi_p2_c2[idx] = (Qmax_p2[cell2] - Qbar->p2[cell2]) / (Qface_c2->p2[idx] - Qbar->p2[cell2]);
        }

        if (Qface_c2->p3[idx] < Qmin_p3[cell2]) {
            phi_p3_c2[idx] = (Qmin_p3[cell2] - Qbar->p3[cell2]) / (Qface_c2->p3[idx] - Qbar->p3[cell2]);
        } else if (Qface_c2->p3[idx] > Qmax_p3[cell2]) {
            phi_p3_c2[idx] = (Qmax_p3[cell2] - Qbar->p3[cell2]) / (Qface_c2->p3[idx] - Qbar->p3[cell2]);
        }

        if (Qface_c2->p4[idx] < Qmin_p4[cell2]) {
            phi_p4_c2[idx] = (Qmin_p4[cell2] - Qbar->p4[cell2]) / (Qface_c2->p4[idx] - Qbar->p4[cell2]);
        } else if (Qface_c2->p4[idx] > Qmax_p4[cell2]) {
            phi_p4_c2[idx] = (Qmax_p4[cell2] - Qbar->p4[cell2]) / (Qface_c2->p4[idx] - Qbar->p4[cell2]);
        }
    }

    // Find the min phi in each cell
    std::vector<double> phi_p1(size->num_cells, 1), phi_p2(size->num_cells, 1), phi_p3(size->num_cells, 1), phi_p4(size->num_cells, 1);

    for (int idx = 0; idx < size->num_faces; idx ++) {
        int cell1 = mesh->face_cell1[idx];
        int cell2 = mesh->face_cell2[idx];

        phi_p1[cell1] = std::min({ phi_p1_c1[idx] , phi_p1[cell1] });
        phi_p4[cell1] = std::min({ phi_p4_c1[idx] , phi_p4[cell1] });
        phi_p2[cell1] = std::min({ phi_p2_c1[idx] , phi_p2[cell1] });
        phi_p3[cell1] = std::min({ phi_p3_c1[idx] , phi_p3[cell1] });

        if (cell2 >= 0) {
            phi_p1[cell2] = std::min({ phi_p1_c2[idx] , phi_p1[cell2] });
            phi_p4[cell2] = std::min({ phi_p4_c2[idx] , phi_p4[cell2] });
            phi_p2[cell2] = std::min({ phi_p2_c2[idx] , phi_p2[cell2] });
            phi_p3[cell2] = std::min({ phi_p3_c2[idx] , phi_p3[cell2] });
        }
    }

    // Limit the gradients 
    for (int idx = 0; idx < size->num_cells; idx ++) {
        Qbar->Qxp1[idx] = Qbar->Qxp1[idx]*phi_p1[idx];
        Qbar->Qxp2[idx] = Qbar->Qxp2[idx]*phi_p2[idx];
        Qbar->Qxp3[idx] = Qbar->Qxp3[idx]*phi_p3[idx];
        Qbar->Qxp4[idx] = Qbar->Qxp4[idx]*phi_p4[idx];

        Qbar->Qyp1[idx] = Qbar->Qyp1[idx]*phi_p1[idx];
        Qbar->Qyp2[idx] = Qbar->Qyp2[idx]*phi_p2[idx];
        Qbar->Qyp3[idx] = Qbar->Qyp3[idx]*phi_p3[idx];
        Qbar->Qyp4[idx] = Qbar->Qyp4[idx]*phi_p4[idx];
    }
}