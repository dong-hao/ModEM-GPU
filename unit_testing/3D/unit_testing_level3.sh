#!/bin/bash
#
# Forward modeling setup was created using the unit_testing_setup.m script in Matlab. Cartesian and
# spherical models were setup to mimic each other, however, there is no meaningful conversion between
# the two. However, the results should be directly comparable.
#
# Comparison with ModEMM done using FWD_3D.m script in Matlab.
#
# Unit testing level 3: runs the inversion and compares the outputs with existing values computed
# using a stable version of the code. Should NOT change substantially if edits are purely technical,
# but details may vary if the algorithm or covariance settings are modified in any way. 

# SET THE PATH
export DYLD_LIBRARY_PATH=".:/opt/intel/lib:/opt/intel/oneapi/mkl/latest/lib:$DYLD_LIBRARY_PATH"

# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - VERSION NAME, 3 - NUMBER OF CORES
#
if [ $# -gt 0 ]; then
    EXEC=$1
else
	EXEC="../../f90/Mod3DMT"
fi

if [ $# -gt 1 ]; then
	version=$2
else
	version='working'
fi

if [ $# -gt 2 ]; then
	ncores=$3
else
	ncores=1
fi

test_dir=$(date '+%Y-%m-%d')-$version
echo "Testing $EXEC in directory $test_dir"

# CREATE TEST OUTPUT FOLDER
mkdir -p OUTPUT
mkdir -p "OUTPUT/$test_dir"
#
# ENTER TEST OUTPUT FOLDER
cd "OUTPUT/$test_dir"
#
#
cmd="../../$EXEC -I NLCG ../../Halfspace+Spherical.rho ../../CylinderModel_6per_Matlab_noisy.dat ../../CylinderModel.inv ../../CylinderModel.fwd"
#
# 
if [ $# -gt 2 ]; then
#
echo "#### START FORWARD MPI TEST WITH $ncores CORES AT $(date)" | tee level3.out
#
#
echo "#### COMMAND LINE: mpirun -n $ncores $cmd "  | tee -a level3.out
#
#
eval mpirun -n $ncores $cmd | tee -a level3.out
#
else
#
echo "#### START FORWARD LOCAL TEST AT $(date)" | tee level3.out
#
#
echo "#### COMMAND LINE: $cmd "  | tee -a level3.out
#
#
eval $cmd | tee -a level3.out
#
fi

#
cd ..
#
#
cd ..
#
#
exit 0
#
# END OF SCRIPT