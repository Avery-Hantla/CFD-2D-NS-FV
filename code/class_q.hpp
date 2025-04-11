#ifndef CLASS_Q_HPP
    #define CLASS_Q_HPP
    #include <vector>
    #include <cmath>
    #include "class_flow.hpp"

    class class_Q {
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

            // Vn avg for entire face!! Only used for q face
            std::vector<double>Vn_avg;
            std::vector<double>c_avg;

            std::vector<double> Qxp1;
            std::vector<double> Qxp2;
            std::vector<double> Qxp3;
            std::vector<double> Qxp4;

            std::vector<double> Qyp1;
            std::vector<double> Qyp2;
            std::vector<double> Qyp3;
            std::vector<double> Qyp4;

            double gamma;

            void init(int size, class_flow* flow) {
                gamma = flow->gamma;
                P.assign(size,flow->P);
                rho.assign(size,flow->rho);
                u.assign(size,flow->u);
                v.assign(size,flow->v);

                Vn_avg.assign(size,std::sqrt(flow->v*flow->v+flow->u*flow->u));
                c_avg.assign(size,std::sqrt(gamma*(flow->P/flow->rho)));

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

                Qxp1.assign(size, 0);
                Qxp2.assign(size, 0);
                Qxp3.assign(size, 0);
                Qxp4.assign(size, 0);

                Qyp1.assign(size, 0);
                Qyp2.assign(size, 0);
                Qyp3.assign(size, 0);
                Qyp4.assign(size, 0);
            }       

            void init(int size, double gamma_in) {
                gamma = gamma_in;
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

                Vn_avg.assign(size,-101);
                c_avg.assign(size,-101);

                p1.assign(size,-101);
                p2.assign(size,-101);
                p3.assign(size,-101);
                p4.assign(size,-101);

                Qxp1.assign(size, 0);
                Qxp2.assign(size, 0);
                Qxp3.assign(size, 0);
                Qxp4.assign(size, 0);

                Qyp1.assign(size, 0);
                Qyp2.assign(size, 0);
                Qyp3.assign(size, 0);
                Qyp4.assign(size, 0);
            }   

            void updateflow() {
                rho = p1;
                E = p4;
                for (int idx = 0; idx < p1.size(); idx++) {
                    u[idx] = p2[idx]/rho[idx];
                    v[idx] = p3[idx]/rho[idx];
                    P[idx] = (E[idx]-(0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx])))*(gamma-1);
                    c[idx] = std::sqrt(gamma*(P[idx]/rho[idx]));
                }
            }

            void updateflow_int(int idx) {
                rho[idx] = p1[idx];
                E[idx] = p4[idx];
                u[idx] = p2[idx]/rho[idx];
                v[idx] = p3[idx]/rho[idx];
                P[idx] = (E[idx]-(0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx])))*(gamma-1);
                c[idx] = std::sqrt(gamma*(P[idx]/rho[idx]));
            }

            void updateQ() {
                p1 = rho;
                for (int idx = 0; idx < p1.size(); idx++) {
                    p2[idx] = rho[idx]*u[idx];
                    p3[idx] = rho[idx]*v[idx];
                    E[idx] = (P[idx]/(gamma-1)) + 0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx]);
                }
                p4 = E;
            }

            int findQf(int cell_num, int face_num) { // Used to find Qi face (NOT CELL Q)
                return cell_num*6 + face_num;
            }
    };
#endif
