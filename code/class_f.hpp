#ifndef CLASS_F_HPP
    #define CLASS_F_HPP
    #include <vector>
    #include "class_mesh.hpp"
    #include "class_q.hpp"
    #include "class_flow.hpp"
    #include "struct_size.hpp"

    class class_F {
        private: 
            double gamma;

        public:
            std::vector<double>p1; 
            std::vector<double>p2;
            std::vector<double>p3; 
            std::vector<double>p4; 

            std::vector<double>Vn;

            void init(int size, double gamma_in) {
                p1.assign(size, -101);
                p2.assign(size, -101);
                p3.assign(size, -101);
                p4.assign(size, -101);
                Vn.assign(size, -101);
                gamma = gamma_in;
            }

            void update(class_mesh* mesh, class_Q* Qface, struct_size* size) { // update ALL
                for (int idx; idx < p1.size(); idx++) {
                    Vn[idx] = Qface->u[idx] * mesh->face_nx[idx] + Qface->v[idx] * mesh->face_ny[idx];

                    p1[idx] = Qface->rho[idx]*Vn[idx];
                    p2[idx] = Qface->rho[idx]*Qface->u[idx]*Vn[idx] + Qface->P[idx]*mesh->face_nx[idx];
                    p3[idx] = Qface->rho[idx]*Qface->v[idx]*Vn[idx] + Qface->P[idx]*mesh->face_ny[idx];
                    p4[idx] = Vn[idx]*(Qface->E[idx] + Qface->P[idx]);
                }
            }
    };
#endif