#include <vector>
#include <cmath>
#include <algorithm>

#include "class_mesh.hpp"
#include "struct_size.hpp"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void minmod(class_mesh* mesh, struct_size* size, std::vector<double> &slope1x, std::vector<double> &slope2x, std::vector<double> &slope3x, std::vector<double> &slope4x, std::vector<double> &slope1y, std::vector<double> &slope2y, std::vector<double> &slope3y, std::vector<double> &slope4y) {
    std::vector<double> minSlopex_p1(size->num_faces, -101), minSlopex_p2(size->num_faces, -101), minSlopex_p3(size->num_faces, -101), minSlopex_p4(size->num_faces, -101);
    std::vector<double> minSlopey_p1(size->num_faces, -101), minSlopey_p2(size->num_faces, -101), minSlopey_p3(size->num_faces, -101), minSlopey_p4(size->num_faces, -101);

    for (int idx = 0; idx < size->num_faces; idx ++) {
        int cell1 = mesh->face_cell1[idx];
        int cell2 = mesh->face_cell2[idx];
        int skipc1 = 0, skipc2 = 0;

        (minSlopex_p1[cell1] == -101) ? minSlopex_p1[idx] = slope1x[idx]: false;
        (minSlopex_p2[cell1] == -101) ? minSlopex_p2[idx] = slope2x[idx]: false;
        (minSlopex_p3[cell1] == -101) ? minSlopex_p3[idx] = slope3x[idx]: false;
        (minSlopex_p4[cell1] == -101) ? minSlopex_p4[idx] = slope4x[idx]: skipc1 = 1;

        (minSlopey_p1[cell1] == -101) ? minSlopey_p1[idx] = slope1y[idx]: false;
        (minSlopey_p2[cell1] == -101) ? minSlopey_p2[idx] = slope2y[idx]: false;
        (minSlopey_p3[cell1] == -101) ? minSlopey_p3[idx] = slope3y[idx]: false;
        (minSlopey_p4[cell1] == -101) ? minSlopey_p4[idx] = slope4y[idx]: skipc1 = 1;

        
        
        if (skipc1 == 0) {
            minSlopex_p1[cell1] = sgn(minSlopex_p1[cell1]) * std::min({ std::abs(minSlopex_p1[cell1]) , std::abs(slope1x[idx]) });
            minSlopex_p2[cell1] = sgn(minSlopex_p2[cell1]) * std::min({ std::abs(minSlopex_p2[cell1]) , std::abs(slope2x[idx]) });
            minSlopex_p3[cell1] = sgn(minSlopex_p3[cell1]) * std::min({ std::abs(minSlopex_p3[cell1]) , std::abs(slope3x[idx]) });
            minSlopex_p4[cell1] = sgn(minSlopex_p4[cell1]) * std::min({ std::abs(minSlopex_p4[cell1]) , std::abs(slope4x[idx]) });

            minSlopey_p1[cell1] = sgn(minSlopey_p1[cell1]) * std::min({ std::abs(minSlopey_p1[cell1]) , std::abs(slope1y[idx]) });
            minSlopey_p2[cell1] = sgn(minSlopey_p2[cell1]) * std::min({ std::abs(minSlopey_p2[cell1]) , std::abs(slope2y[idx]) });
            minSlopey_p3[cell1] = sgn(minSlopey_p3[cell1]) * std::min({ std::abs(minSlopey_p3[cell1]) , std::abs(slope3y[idx]) });
            minSlopey_p4[cell1] = sgn(minSlopey_p4[cell1]) * std::min({ std::abs(minSlopey_p4[cell1]) , std::abs(slope4y[idx]) });
        }

        if (cell2 >= 0) {
            (minSlopex_p1[cell2] == -101) ? minSlopex_p1[idx] = slope1x[idx]: false;
            (minSlopex_p2[cell2] == -101) ? minSlopex_p2[idx] = slope2x[idx]: false;
            (minSlopex_p3[cell2] == -101) ? minSlopex_p3[idx] = slope3x[idx]: false;
            (minSlopex_p4[cell2] == -101) ? minSlopex_p4[idx] = slope4x[idx]: skipc2 = 1;

            (minSlopey_p1[cell2] == -101) ? minSlopey_p1[idx] = slope1y[idx]: false;
            (minSlopey_p2[cell2] == -101) ? minSlopey_p2[idx] = slope2y[idx]: false;
            (minSlopey_p3[cell2] == -101) ? minSlopey_p3[idx] = slope3y[idx]: false;
            (minSlopey_p4[cell2] == -101) ? minSlopey_p4[idx] = slope4y[idx]: skipc2 = 1;

            if (skipc2 == 0) {
                minSlopex_p1[cell2] = sgn(minSlopex_p1[cell2]) * std::min({ std::abs(minSlopex_p1[cell2]) , std::abs(slope1x[idx]) });
                minSlopex_p2[cell2] = sgn(minSlopex_p2[cell2]) * std::min({ std::abs(minSlopex_p2[cell2]) , std::abs(slope2x[idx]) });
                minSlopex_p3[cell2] = sgn(minSlopex_p3[cell2]) * std::min({ std::abs(minSlopex_p3[cell2]) , std::abs(slope3x[idx]) });
                minSlopex_p4[cell2] = sgn(minSlopex_p4[cell2]) * std::min({ std::abs(minSlopex_p4[cell2]) , std::abs(slope4x[idx]) });

                minSlopey_p1[cell2] = sgn(minSlopey_p1[cell2]) * std::min({ std::abs(minSlopey_p1[cell2]) , std::abs(slope1y[idx]) });
                minSlopey_p2[cell2] = sgn(minSlopey_p2[cell2]) * std::min({ std::abs(minSlopey_p2[cell2]) , std::abs(slope2y[idx]) });
                minSlopey_p3[cell2] = sgn(minSlopey_p3[cell2]) * std::min({ std::abs(minSlopey_p3[cell2]) , std::abs(slope3y[idx]) });
                minSlopey_p4[cell2] = sgn(minSlopey_p4[cell2]) * std::min({ std::abs(minSlopey_p4[cell2]) , std::abs(slope4y[idx]) });
            }
        }
    }

    // Put minmum cell slope in faces
    for (int idx = 0; idx < size->num_cells; idx ++) {
        for (int jdx = 0; jdx < mesh->cell_face_count[idx]; jdx ++) {
            int face_num = mesh->find_cell_face(idx,jdx);
            slope1x[face_num] = minSlopex_p1[idx];
            slope2x[face_num] = minSlopex_p2[idx];
            slope3x[face_num] = minSlopex_p3[idx];
            slope4x[face_num] = minSlopex_p4[idx];

            slope1y[face_num] = minSlopey_p1[idx];
            slope2y[face_num] = minSlopey_p2[idx];
            slope3y[face_num] = minSlopey_p3[idx];
            slope4y[face_num] = minSlopey_p4[idx];
        }
    }

}