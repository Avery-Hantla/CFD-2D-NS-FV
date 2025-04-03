#ifndef CLASS_FLOW_HPP
  #define CLASS_FLOW_HPP
  class class_flow{
    public:
      double gamma;
      double P; 
      double rho;
      double u;
      double v;
      double E;

      double Q1;
      double Q2;
      double Q3;
      double Q4;

      void updateQ() {
        Q1 = rho;
        Q2 = rho*u;
        Q3 = rho*v;
        E = (P/(gamma-1)) + 0.5*rho*(u*u+v*v);
        Q4 = E;
      }
  };
#endif