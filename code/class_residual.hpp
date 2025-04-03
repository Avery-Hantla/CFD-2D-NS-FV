#ifndef CLASS_RESIDUAL_HPP
    #define CLASS_RESIDUAL_HPP
    #include <vector>

    class class_residual{
        public:
            std::vector<double> p1;
            std::vector<double> p2;
            std::vector<double> p3;
            std::vector<double> p4;

            void init(int size) {
                p1.assign(size, 0);
                p2.assign(size, 0);
                p3.assign(size, 0);
                p4.assign(size, 0);
            }

            void reset() {
                std::fill(p1.begin(), p1.end(), 0);
                std::fill(p2.begin(), p2.end(), 0);
                std::fill(p3.begin(), p3.end(), 0);
                std::fill(p4.begin(), p4.end(), 0);
            }
    };
#endif