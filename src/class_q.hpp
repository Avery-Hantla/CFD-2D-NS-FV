#ifndef CLASS_Q_HPP
    #define CLASS_Q_HPP
    #include <vector>
    #include <cmath>
    #include "class_flow.hpp"

    class class_Q {
        public:
            double R;
            double gamma;
            double k;
            double mu;

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

            std::vector<double>T;  

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

            // Taus for viscous flow on faces
            std::vector<double> tau_xx;
            std::vector<double> tau_xy; 
            std::vector<double> tau_yy; 

            std::vector<double> dTdx;
            std::vector<double> dTdy;

            void init(int size, class_flow* flow) {
                gamma = flow->gamma;
                R = flow->R;
                mu = flow->mu;
                k = flow->k;

                P.assign(size,flow->P);
                rho.assign(size,flow->rho);
                u.assign(size,flow->u);
                v.assign(size,flow->v);

                Vn_avg.assign(size,std::sqrt(flow->v*flow->v+flow->u*flow->u));
                c_avg.assign(size,std::sqrt(gamma*(flow->P/flow->rho)));

                E.assign(size,-101);
                c.assign(size,-101);
                T.assign(size,-101);
                for (int idx = 0; idx < size; idx ++) {
                    E[idx] = (P[idx]/(gamma-1)) + 0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx]);
                    c[idx] = std::sqrt(gamma*(P[idx]/rho[idx]));
                    T[idx] = P[idx]/(rho[idx]*R);
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

                tau_xx.assign(size, 0);
                tau_xy.assign(size, 0);
                tau_yy.assign(size, 0);

                dTdx.assign(size, 0);
                dTdy.assign(size, 0);
            }       

            void init(int size, double gamma_in, double R_in) {
                gamma = gamma_in;
                R = R_in;
                P.assign(size,-101);
                rho.assign(size,-101);
                u.assign(size,-101);
                v.assign(size,-101);
                E.assign(size,-101);
                c.assign(size,-101);
                T.assign(size,-101);

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
                    T[idx] = P[idx]/(rho[idx]*R);
                }
            }

            void updateflow_int(int idx) {
                rho[idx] = p1[idx];
                E[idx] = p4[idx];
                u[idx] = p2[idx]/rho[idx];
                v[idx] = p3[idx]/rho[idx];
                P[idx] = (E[idx]-(0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx])))*(gamma-1);
                c[idx] = std::sqrt(gamma*(P[idx]/rho[idx]));
                T[idx] = P[idx]/(rho[idx]*R); 
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

            void comp_vis(int idx, double ux, double uy, double vx, double vy, double Px, double Py, double rhox, double rhoy) {
                    tau_xx[idx] = 2*mu*(ux - (ux+vy)/3);
                    tau_xy[idx] = mu*(uy+vx);
                    tau_yy[idx] = 2*mu*(vy - (ux+vy)/3);
                    dTdx[idx] = (Px*rho[idx] - rhox*P[idx])/(R*rho[idx]*rho[idx]);
                    dTdy[idx] = (Py*rho[idx] - rhoy*P[idx])/(R*rho[idx]*rho[idx]);
            }

            // int findQf(int cell_num, int face_num) { // Used to find Qi face (NOT CELL Q)
            //     return cell_num*6 + face_num;
            // }
    };
#endif
