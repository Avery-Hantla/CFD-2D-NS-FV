#ifndef CLASS_Q_HPP
    #define CLASS_Q_HPP
    #include <vector>
    #include <cmath>
    #include "class_flow.hpp"

    class class_Q {
        private:
            double gamma;

        public:
            std::vector<double>p1; 
            std::vector<double>p2;
            std::vector<double>p3; 
            std::vector<double>p4;  

            std::vector<double>rho;
            std::vector<double>u;
            std::vector<double>v;
            std::vector<double>E;
            std::vector<double>P;
            std::vector<double>c;  

            void init(int size, class_flow* flow) {
                gamma = flow->gamma;
                P.assign(size,flow->P);
                rho.assign(size,flow->rho);
                u.assign(size,flow->u);
                v.assign(size,flow->v);

                E.assign(size,-101);
                c.assign(size,-101);
                for (int idx = 0; idx < size; idx ++) {
                    E[idx] = (P[idx]/(gamma-1)) + 0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx]);
                    c[idx] = std::sqrt(gamma*(P[idx]/rho[idx]));
                }
                
                p1.assign(size,-101);
                p2.assign(size,-101);
                p3.assign(size,-101);
                p4.assign(size,-101);
            }       

            void init(int size) {
                gamma = -101;
                P.assign(size,-101);
                rho.assign(size,-101);
                u.assign(size,-101);
                v.assign(size,-101);
                E.assign(size,-101);
                c.assign(size,-101);

                p1.assign(size,-101);
                p2.assign(size,-101);
                p3.assign(size,-101);
                p4.assign(size,-101);
            }   

            void updateflow() {
                rho = p1;
                E = p4;
                for (int idx; idx < p1.size(); idx++) {
                    u[idx] = p2[idx]/rho[idx];
                    v[idx] = p3[idx]/rho[idx];
                    P[idx] = (E[idx]-0.5*rho[idx]*(u[idx]+v[idx]))*(gamma-1);
                    c[idx] = std::sqrt(gamma*(P[idx]/rho[idx]));
                }
            }

            void updateQ() {
                p1 = rho;
                p4 = E;
                for (int idx; idx < p1.size(); idx++) {
                    p2[idx] = rho[idx]*u[idx];
                    p3[idx] = rho[idx]*v[idx];
                }
            }

            int findQf(int cell_num, int face_num) { // Used to find Qi face (NOT CELL Q)
                return cell_num*6 + face_num;
            }
    };
#endif
