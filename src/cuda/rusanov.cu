/////////////////////////////////////////////////////////
//           Function To Compute Rusanov Flux
/////////////////////////////////////////////////////////

#include "../support/class_q.hpp"
#include "../support/class_mesh.hpp"
#include "../support/class_f.hpp"
#include "../support/class_flow.hpp"
#include "../support/struct_size.hpp"
#include "../support/struct_BC.hpp"

// #include "../solver/rusanov.hpp"

#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>
// #include "cuda.cuh"

// temp move to function in main.cpp
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

        // printf("Hello from block: %u, thread: %u\n", blockIdx.x, threadIdx.x);
    }
}

void rusanov(class_F* F_rusanov, class_Q* Qface_c1, class_Q* Qface_c2, class_mesh* mesh, class_flow* freestream, struct_size* size, struct_BC* BC) {
    class_F F_c1, F_c2;
    F_c1.init(size->num_faces, freestream->gamma);
    F_c2.init(size->num_faces, freestream->gamma);
    
    F_c1.update(mesh, Qface_c1, size); 
    F_c2.update(mesh, Qface_c2, size);

    // Device pointers for vectors
    double *d_cl_Vn, *d_c2_Vn, *d_c1_c, *d_c2_c;
    double *d_Vn_avg, *d_c_avg;

    double *d_Q_c1_p1, *d_Q_c1_p2, *d_Q_c1_p3, *d_Q_c1_p4; 
    double *d_Q_c2_p1, *d_Q_c2_p2, *d_Q_c2_p3, *d_Q_c2_p4;

    double *d_F_c1_p1, *d_F_c1_p2, *d_F_c1_p3, *d_F_c1_p4;
    double *d_F_c2_p1, *d_F_c2_p2, *d_F_c2_p3, *d_F_c2_p4;

    double *d_F_rusanov_p1, *d_F_rusanov_p2, *d_F_rusanov_p3, *d_F_rusanov_p4;

    // Allocate memory on the device for vectors
    cudaMalloc((void **)&d_cl_Vn, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_c2_Vn, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_c1_c, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_c2_c, size->num_faces * sizeof(double));

    cudaMalloc((void **)&d_Vn_avg, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_c_avg, size->num_faces * sizeof(double));
    
    cudaMalloc((void **)&d_Q_c1_p1, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_Q_c1_p2, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_Q_c1_p3, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_Q_c1_p4, size->num_faces * sizeof(double));

    cudaMalloc((void **)&d_Q_c2_p1, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_Q_c2_p2, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_Q_c2_p3, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_Q_c2_p4, size->num_faces * sizeof(double));

    cudaMalloc((void **)&d_F_c1_p1, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_c1_p2, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_c1_p3, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_c1_p4, size->num_faces * sizeof(double));

    cudaMalloc((void **)&d_F_c2_p1, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_c2_p2, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_c2_p3, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_c2_p4, size->num_faces * sizeof(double));

    cudaMalloc((void **)&d_F_rusanov_p1, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_rusanov_p2, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_rusanov_p3, size->num_faces * sizeof(double));
    cudaMalloc((void **)&d_F_rusanov_p4, size->num_faces * sizeof(double));

    // Copy vectors from host to device
    cudaMemcpy(d_cl_Vn, &F_c1.Vn[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c2_Vn, &F_c2.Vn[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c1_c, &Qface_c1->c[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c2_c, &Qface_c2->c[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_Q_c1_p1, &Qface_c1->p1[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q_c1_p2, &Qface_c1->p2[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q_c1_p3, &Qface_c1->p3[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q_c1_p4, &Qface_c1->p4[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_Q_c2_p1, &Qface_c2->p1[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q_c2_p2, &Qface_c2->p2[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q_c2_p3, &Qface_c2->p3[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q_c2_p4, &Qface_c2->p4[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_F_c1_p1, &F_c1.p1[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_c1_p2, &F_c1.p2[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_c1_p3, &F_c1.p3[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_c1_p4, &F_c1.p4[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_F_c2_p1, &F_c2.p1[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_c2_p2, &F_c2.p2[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_c2_p3, &F_c2.p3[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_c2_p4, &F_c2.p4[0], size->num_faces * sizeof(double), cudaMemcpyHostToDevice);

    int nThreads = maxThreadsperblock();
    int nBlocks = ceil((double)size->num_faces / nThreads);

    // Run Kernel
    d_rusanov <<<nBlocks, nThreads>>> (size->num_faces, d_cl_Vn, d_c2_Vn, d_c1_c, d_c2_c, d_Vn_avg, d_c_avg, d_Q_c1_p1, d_Q_c1_p2, d_Q_c1_p3, d_Q_c1_p4, d_Q_c2_p1, d_Q_c2_p2, d_Q_c2_p3, d_Q_c2_p4, d_F_c1_p1, d_F_c1_p2, d_F_c1_p3, d_F_c1_p4, d_F_c2_p1, d_F_c2_p2, d_F_c2_p3, d_F_c2_p4, d_F_rusanov_p1, d_F_rusanov_p2, d_F_rusanov_p3, d_F_rusanov_p4);

    // Check for error
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }

    cudaDeviceSynchronize();

    // Copy vectors from device to host
    cudaMemcpy(&F_rusanov->p1[0], d_F_rusanov_p1, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&F_rusanov->p2[0], d_F_rusanov_p2, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&F_rusanov->p3[0], d_F_rusanov_p3, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&F_rusanov->p4[0], d_F_rusanov_p4, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);

    cudaMemcpy(&Qface_c1->Vn_avg[0], d_Vn_avg, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&Qface_c2->Vn_avg[0], d_Vn_avg, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&Qface_c1->c_avg[0], d_c_avg, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&Qface_c2->c_avg[0], d_c_avg, size->num_faces * sizeof(float), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaDeviceReset();
}
