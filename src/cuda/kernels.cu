#ifndef KERNELS_CU
    #define KERNELS_CU

    #include <cuda_runtime.h>
    #include <iostream>

    __global__ void d_rusanov(int N, double d_cl_Vn[], double d_c2_Vn[], double d_c1_c[], double d_c2_c[], double d_Vn_avg[], double d_c_avg[], double d_Q_c1_p1[], double d_Q_c1_p2[], double d_Q_c1_p3[], double d_Q_c1_p4[], double d_Q_c2_p1[], double d_Q_c2_p2[], double d_Q_c2_p3[], double d_Q_c2_p4[], double d_F_c1_p1[], double d_F_c1_p2[], double d_F_c1_p3[], double d_F_c1_p4[], double d_F_c2_p1[], double d_F_c2_p2[], double d_F_c2_p3[], double d_F_c2_p4[], double d_F_rusanov_p1[], double d_F_rusanov_p2[], double d_F_rusanov_p3[], double d_F_rusanov_p4[]) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        if (idx < N) {
            // Find Vn and c averag 
            double Vn_avg = (d_cl_Vn[idx] + d_c2_Vn[idx])/2;
            double c_avg = (d_c1_c[idx] + d_c2_c[idx])/2;

            // Save for time step calculation
            d_Vn_avg[idx] = Vn_avg;
            d_c_avg[idx] = c_avg;

            // Compute Rusanov Flux
            d_F_rusanov_p1[idx] = (d_F_c1_p1[idx]+d_F_c2_p1[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (d_Q_c2_p1[idx] - d_Q_c1_p1[idx]);
            d_F_rusanov_p2[idx] = (d_F_c1_p2[idx]+d_F_c2_p2[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (d_Q_c2_p2[idx] - d_Q_c1_p2[idx]);
            d_F_rusanov_p3[idx] = (d_F_c1_p3[idx]+d_F_c2_p3[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (d_Q_c2_p3[idx] - d_Q_c1_p3[idx]);
            d_F_rusanov_p4[idx] = (d_F_c1_p4[idx]+d_F_c2_p4[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (d_Q_c2_p4[idx] - d_Q_c1_p4[idx]);

            printf("Hello from block: %u, thread: %u\n", blockIdx.x, threadIdx.x);
        }
    }

#endif