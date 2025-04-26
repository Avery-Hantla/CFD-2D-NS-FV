#ifndef CUDA_HPP
    #define CUDA_HPP
    
    #include <cuda_runtime.h>

    // #include "cuda/kernels.cu"
    // __global__ void d_rusanov(int N, double d_cl_Vn[], double d_c2_Vn[], double d_c1_c[], double d_c2_c[], double d_Vn_avg[], double d_c_avg[], double d_Q_c1_p1[], double d_Q_c1_p2[], double d_Q_c1_p3[], double d_Q_c1_p4[], double d_Q_c2_p1[], double d_Q_c2_p2[], double d_Q_c2_p3[], double d_Q_c2_p4[], double d_F_c1_p1[], double d_F_c1_p2[], double d_F_c1_p3[], double d_F_c1_p4[], double d_F_c2_p1[], double d_F_c2_p2[], double d_F_c2_p3[], double d_F_c2_p4[], double d_F_rusanov_p1[], double d_F_rusanov_p2[], double d_F_rusanov_p3[], double d_F_rusanov_p4[]);

    int maxThreadsperblock() {
        // current CUDA device ID
        int device;

        // object to store properties of the CUDA device
        cudaDeviceProp props;

        // ID of the currently active CUDA device
        cudaGetDevice(&device);

        // Retrieve properties and store in 'props'
        cudaGetDeviceProperties(&props, device);

        // Return max number of threads per block
        return props.maxThreadsPerBlock;
    }
#endif