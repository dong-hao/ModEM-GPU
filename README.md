# ModEM-GPU
Snapshot Repository for the Manuscript "Hybrid CPU-GPU solution to regularized divergence-free curl-curl equations for electromagnetic inversion problems" submitted to CnG

This is a new experimental branch of ModEM that utilize a hybrid CPU-GPU algorithm to rapidly compute the forward and adjoint systems in EM inversion problems. Please see the manuscript for the detailed mathematical formulas and experiments that demonstrates the efficiency of the new method.

## Below is the commit log from the original ModEM SVN repo

introducing the new GPU solvers – 

this (hopefully) will not affect the behavior of CPU solvers. for now, I only implemented the GPU solver in the SP2 version - which is the most efficient CPU version so far. 

Depending on your hardware, the GPU vs CPU speed-up can be anywhere between a few times (for old/weak GPUs) to tens of times for professional ones. 

The code has been tested on various computers from my old laptop with a consumer GPU (bought in 2016) to a brand-new multi-GPU workstation, with GCC+Gfortran 7/9 and CUDA 10/11

Also note that this apparently only works for NVIDIA cards, as the
implementation is through CUDA. Other "universal" interfaces like Kokkos or OpenCL may be prefered for other GPUs - however, those interfaces do not provide the ability to implement kernel-level improvements, yet. 

Those who want to use the code are welcomed to try - although one should bear in mind that the GPU libs and CUDA are rather user-hostile to set-up. You have been warned... See Makefile.gpu file for some basic examples on how to compile the code.

Hao

## TODO
There is still no means to automatically configure the GPU version with the f90/CONFIG/ scripts, as the code now involves a CUDA cpp compiler (NVCC) - need to do it by hand for now - i.e. one needs to set up the $CC and $C_FLAG thingy by corresponding values from a particular system. (need to investigate a new way to do the auto-configuration)

NOTE: this repo will be moved to the USGS code website, once Anna finishes setting up the license and repository there. 

## CONTACT
For more information about this code, or if you have any queries or advices, please do not hesitate to contact me (or Gary/Anna/Naser). 

Email: donghao@cugb.edu.cn; 
       gary.egbert@oregonstate.edu;
       akelbert@usgs.gov;
       meqbel@modem-geophysics.com
