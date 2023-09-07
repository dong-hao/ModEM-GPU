#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <nvml.h>
#include <complex.h>

// simple kernel function that converts double vectors to single
__global__ void real64to32(const double *in, float *out, const int N)
{
    // a position every 64 bits
    // int pos = blockDim.x * blockIdx.x + threadIdx.x;
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	float sing = __double2float_rn(in[pos]);
	out[pos] = sing;
    }
}

// simple kernel function that converts single vectors to double
__global__ void real32to64(const uint32_t *in, uint64_t *out, const int N)
{
    // a position every 32 bits
    // int pos = blockDim.x * blockIdx.x + threadIdx.x;
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	// last sign bit (1)
	uint64_t s = in[pos] & 0x80000000;
	// exponent bits (8)
	uint64_t e = ((in[pos] & 0x7f800000) >> 23);
	e = e + 896;
	// mantissa bits (23)
	uint64_t m = in[pos] & 0x7fffff;
	// double through bitwise or
	uint64_t doub = (s<<32) | (e<<52) | (m<<29);
        // a new position every 64 bits
	out[pos] = doub;
    }
 }

__global__ void hada_real(const double *ina, const double *inb, double *out, const int N)
{
    //hadamard multiplication, as real 
    // c = a * b
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	out[pos] = ina[pos]*inb[pos];
    }
 }

__global__ void hada_cmplx(const double *ina, const double *inb, double *out, const int N)
{
    //hadamard multiplication, as complex
    //(a + bi) * (c + di) = (ac -bd) + (ad + bc)i
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
        if ( pos%2 == 0) // 0 2 4 
        {
	        out[pos] = ina[pos]*inb[pos]-ina[pos+1]*inb[pos+1];
        }
        else // 1 3 5 
        {
	        out[pos] = ina[pos-1]*inb[pos]+ina[pos]*inb[pos-1];
        }
    }
 }

__global__ void xpby_real(const double *x, const double b, double *y, const int N)
{
    // CUDA kernel implementing xpby:
    // y = x + b*y
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	        y[pos] = x[pos]+b*y[pos];
    }
}

__global__ void xpby_cmplx(const double *x, double br, double bi, double *y, const int N)
{
    // CUDA kernel implementing xpby:
    // y = x + b*y
    // complex makes this a little tricky - need a temp variable
    int pos = ((blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x)*2;
    if (pos<N) 
    {  
        // lower 64bit --> real
	    double lower = x[pos] + br*y[pos] - bi*y[pos+1];
        // higher 64bit --> imag
	    y[pos+1] = x[pos+1] + br*y[pos+1] + bi*y[pos];
        // and copy them back
        y[pos] = lower;
    }
}

__global__ void p_update_real(const double *v, const double *r, const double beta, const double omega, double *p, const int N)
{
    // CUDA kernel implementing p update:
    // p = r + beta * (p - omega*v)
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	    p[pos] = r[pos] + beta*(p[pos] - omega*v[pos]);
    }
}

__global__ void p_update_cmplx(const double *r, const double *v, const double br, const double bi, const double wr, const double wi, double *p, const int N)
{
    // CUDA kernel implementing p update:
    // p = r + beta * (p - omega*v)
    int pos = ((blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x)*2;
    if (pos<N) 
    {  
        // lower 64bit --> real p - omega*v
	    double lower = p[pos] - wr*v[pos] + wi*v[pos+1];
        // higher 64bit --> imag p - omega*v
	    double higher = p[pos+1] - wr*v[pos+1] - wi*v[pos];
        // lower 64bit --> real r + beta * (p - omega*v)
	    p[pos] = r[pos] + br*lower - bi*higher;
        // higher 64bit --> imag r + beta * (p - omega*v)
	    p[pos+1] = r[pos+1] + br*higher + bi*lower;
    }
}

__global__ void x_update_cmplx(const double *p, const double *s, const double ar, const double ai, const double wr, const double wi, double *x, const int N)
{
    // CUDA kernel implementing p update:
    int pos = ((blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x)*2;
    if (pos<N) 
    {  
        // x = x + alpha * ph + omega * sh
        // lower 64bit --> real x + alpha*ph
	    x[pos] = x[pos] + ar*p[pos] - ai*p[pos+1] + wr*s[pos] - wi*s[pos];
        // higher 64bit --> imag x + alpha*ph
	    x[pos+1]= x[pos+1] + ar*p[pos+1] + ai*p[pos] + wr*s[pos+1] + wi*s[pos];
    }
}

__global__ void reduce_real(double *a, double s, const int N)
{
    // CUDA kernel implementing reduce s = sum(a)
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    int number_of_threads = (blockDim.x*blockDim.y);
    int step = 1; //initial step
    while (number_of_threads > 0)
    {
        if (pos < number_of_threads) 
        {
            const int left = pos * step*2; 
            const int right = left + step;
            a[left] += a[right];
        }
        step <<= 1;
        number_of_threads >>= 1;
        __syncthreads();
    }
}

// function to wait for sometime
void sleep(int seconds)
{
    // Converting time into milli_seconds
    int milli_seconds = 1000 * seconds;
    // Storing start time
    clock_t start_time = clock();
    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds)
        ;
}

// function called from main fortran program
// single to double conversion
extern "C" void kernelc_s2d(const uint32_t *a_d, uint64_t *b_d, int Np)
{
    //uint32_t  *a_d;  // declare GPU vector double (but stored as uint)
    //uint64_t  *b_d;  // declare GPU vector float
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np*2;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    real32to64<<< grids, blocks >>>( a_d, b_d, N);

    return;
}

// function called from main fortran program
// double to single conversion
extern "C" void kernelc_d2s(const double *a_d, float *b_d, int Np)
{
    //double  *a_d;  // declare GPU vector double 
    //float  *b_d;  // declare GPU vector float
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np*2;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    real64to32<<< grids, blocks >>>( a_d, b_d, N);

    return;
}

// function called from main fortran program
// hadamard multiply real version
extern "C" void kernelc_hadar(const double *a_d, const double *b_d, double *c_d, int Np)
{
    //double  *a_d;  // declare GPU vector double 
    //double  *b_d;  // declare GPU vector double
    //double  *c_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    hada_real<<< grids, blocks >>>( a_d, b_d, c_d, N);

    return;
}

// function called from main fortran program
// hadamard multiply complex version
extern "C" void kernelc_hadac(const double *a_d, const double *b_d, double *c_d, int Np)
{
    //double  *a_d;  // declare GPU vector double 
    //double  *b_d;  // declare GPU vector double
    //double  *c_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np*2;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    hada_cmplx<<< grids, blocks >>>( a_d, b_d, c_d, N);

    return;
}

// function called from main fortran program
// xpby complex version
extern "C" void kernelc_xpbyc(const double *x_d, double _Complex b, double *y_d, int Np)
{
    //double  *x_d;  // declare GPU vector double 
    //double  *y_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = 2*Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // get upper 64 bit
    double br = creal(b);
    // get lower 64 bit
    double bi = cimag(b);
    // call function on GPU
    xpby_cmplx<<< grids, blocks >>>( x_d, br, bi, y_d, N);

    return;
}

// function called from main fortran program
// update p complex version
extern "C" void kernelc_update_pc(const double *r_d, const double *v_d, double _Complex beta, double _Complex omega, double *p_d, int Np)
{
    //double  *a_d;  // declare GPU vector double 
    //double  *b_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = 2*Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // get upper 64 bit
    double br = creal(beta);
    // get lower 64 bit
    double bi = cimag(beta);
    // get upper 64 bit
    double wr = creal(omega);
    // get lower 64 bit
    double wi = cimag(omega);
    // call function on GPU
    p_update_cmplx<<< grids, blocks >>>( r_d, v_d, br, bi, wr, wi, p_d, N);

    return;
}

// function called from main fortran program
// update x complex version
extern "C" void kernelc_update_xc(const double *ph_d, const double *sh_d, double _Complex alpha, double _Complex omega, double *x_d, int Np)
{
    //double  *ph_d;  // declare GPU vector double 
    //double  *sh_d;  // declare GPU vector double
    //double  *x_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = 2*Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // get upper 64 bit
    double ar = creal(alpha);
    // get lower 64 bit
    double ai = cimag(alpha);
    // get upper 64 bit
    double wr = creal(omega);
    // get lower 64 bit
    double wi = cimag(omega);
    // call function on GPU
    x_update_cmplx<<< grids, blocks >>>( ph_d, sh_d, ar, ai, wr, wi, x_d, N);

    return;
}

// function called from main fortran program 
extern "C" int kernelc_hookCtx(int dev_idx)
{
    nvmlReturn_t result;
    cudaError_t err;
    CUcontext thisCtx;
    nvmlDevice_t device;
    nvmlPciInfo_t pci_bus;
    nvmlUtilization_t usage;
    unsigned int Ndevice;
    // unsigned int i;
    // unsigned int nProc = 32;
    // nvmlProcessInfo_t pInfo[nProc];
    char device_name[NVML_DEVICE_NAME_BUFFER_SIZE];

    // First initialize NVML library
    result = nvmlInit();
    if (NVML_SUCCESS != result)
    { 
        printf("Failed to initialize NVML: %s\n", nvmlErrorString(result));
        return 1;
    }
    result = nvmlDeviceGetCount(&Ndevice);
    if (NVML_SUCCESS != result)
    { 
        printf("Failed to initialize NVML: %s\n", nvmlErrorString(result));
        return 1;
    }
    if (dev_idx + 1 > Ndevice) 
    { 
        printf("Error: not enough devices detected \n");
        printf("required device idx = %i, Ndevice = %u \n", dev_idx, Ndevice);
        goto Error;
    }
    err = cudaSetDevice(dev_idx);
	if( err != cudaSuccess )
    {
        printf("Failed to set device %i: %u \n", dev_idx, err);
        goto Error;
    }
    while(true)
    {
        // Query for device handle to perform operations on a device
        // You can also query device handle by other features like:
        // nvmlDeviceGetHandleBySerial
        // nvmlDeviceGetHandleByPciBusId
        result = nvmlDeviceGetHandleByIndex(dev_idx, &device);
        if (NVML_SUCCESS != result)
        { 
            printf("Failed to get handle for device %i: %s\n", dev_idx, nvmlErrorString(result));
            goto Error;
        }
        // now get the current context (if any)
        cuCtxGetCurrent(&thisCtx);
        if(thisCtx == NULL) // first call to this device
        {
            cuCtxCreate(&thisCtx, 0, dev_idx);
            result = nvmlDeviceGetName(device, device_name, NVML_DEVICE_NAME_BUFFER_SIZE);
            if (NVML_SUCCESS != result)
            { 
                printf("Failed to get name of device %i: %s\n", dev_idx, nvmlErrorString(result));
                goto Error;
            }
            // pci.busId is very useful to know which device physically 
            // you're talking to
            // Using PCI identifier you can also match nvmlDevice 
            // handle to CUDA device.
            result = nvmlDeviceGetPciInfo(device, &pci_bus);
            if (NVML_SUCCESS != result)
            { 
                printf("Failed to get pci info for device %i: %s\n", dev_idx, nvmlErrorString(result));
                goto Error;
            }
            break;// just go ahead to initialize the case
        }
        else// see if this device is available
        {
            cuCtxSetCurrent(thisCtx);
            result = nvmlDeviceGetUtilizationRates( device, &usage );
            if (NVML_SUCCESS != result)
            {
                printf("Failed to get usage for device %i: %s\n", dev_idx, nvmlErrorString(result));
                goto Error;
            }
            if (usage.memory+usage.gpu < 15)
            {
                printf(" # Dev Status  : GPU-util=%u GPU-mem=%u \n", usage.gpu, usage.memory );
                result = nvmlDeviceGetName(device, device_name, NVML_DEVICE_NAME_BUFFER_SIZE);
                if (NVML_SUCCESS != result)
                { 
                    printf("Failed to get name of device %i: %s\n", dev_idx, nvmlErrorString(result));
                    goto Error;
                }
                // pci.busId is very useful to know which device physically 
                // you're talking to
                // Using PCI identifier you can also match nvmlDevice 
                // handle to CUDA device.
                result = nvmlDeviceGetPciInfo(device, &pci_bus);
                if (NVML_SUCCESS != result)
                { 
                    printf("Failed to get pci info for device %i: %s\n", dev_idx, nvmlErrorString(result));
                    goto Error;
                }
                break;
            }
        }
        sleep(0.05);
    }
    printf(" # Dev Selected:  %i. %s [%s]\n", dev_idx, device_name, pci_bus.busId);
    result = nvmlShutdown();
    if (NVML_SUCCESS != result)
        printf("Failed to shutdown NVML: %s\n", nvmlErrorString(result));

    return 0;
Error:
    result = nvmlShutdown();
    if (NVML_SUCCESS != result)
        printf("Failed to shutdown NVML: %s\n", nvmlErrorString(result));

    return 1;
}

// externel function to tell the number of GPU devices
extern "C" int kernelc_getDevNum()
{
    nvmlReturn_t result;
    unsigned int nGPU;
    int          Ndevice;

    // First initialize NVML library
    result = nvmlInit();
    if (NVML_SUCCESS != result)
    { 
        printf("Failed to initialize NVML: %s\n", nvmlErrorString(result));
        return -1;
    }

    // number of devices
    result = nvmlDeviceGetCount(&nGPU);
    if (NVML_SUCCESS != result)
    { 
        printf("Failed to query device count: %s\n", nvmlErrorString(result));
        goto Error;
    }
    // go silent
    // printf("Found %u device%s\n", Ndevice, Ndevice != 1 ? "s" : " ");

    result = nvmlShutdown();
    if (NVML_SUCCESS != result)
        printf("Failed to shutdown NVML: %s\n", nvmlErrorString(result));

    Ndevice = nGPU;
    return Ndevice;
Error:
    result = nvmlShutdown();
    if (NVML_SUCCESS != result)
        printf("Failed to shutdown NVML: %s\n", nvmlErrorString(result));

    Ndevice = 0;
    return -1;
}
