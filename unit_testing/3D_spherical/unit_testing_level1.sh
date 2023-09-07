#!/bin/bash
#
# Forward modeling setup was created using the unit_testing_setup.m script in Matlab. Cartesian and
# spherical models were setup to mimic each other, however, there is no meaningful conversion between
# the two. However, the results should be directly comparable.
#
# Comparison with ModEMM done using FWD_3D_spherical.m script in Matlab.
#
# Unit testing level 1: runs the forward solver and compares the outputs with existing values computed
# using a stable version of the code. Should NOT change if edits are purely technical. 

# SET THE PATH
export DYLD_LIBRARY_PATH=".:/opt/intel/lib:/opt/intel/oneapi/mkl/latest/lib:$DYLD_LIBRARY_PATH"

# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - VERSION NAME, 3 - NUMBER OF CORES
#
if [ $# -gt 0 ]; then
    EXEC=$1
else
	EXEC="../../f90/Mod3DMTs"
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
cd "OUTPUT"
#
# 
if [ $# -gt 2 ]; then
#
cmd="../$EXEC -F ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_6per_Spherical_Template.dat $test_dir/CylinderModel_6per.dat $test_dir/CylinderModel_6per.soln ../CylinderModel.fwd"
#
#
echo "#### START FORWARD MPI TEST WITH $ncores CORES AT $(date)" | tee $test_dir/level1.out
#
#
echo "#### COMMAND LINE: mpirun -n $ncores $cmd "  | tee -a $test_dir/level1.out
#
#
eval mpirun -n $ncores $cmd | tee -a $test_dir/level1.out
#
else
#
cmd="../$EXEC -F ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_1per_Spherical_Template.dat $test_dir/CylinderModel_1per.dat $test_dir/CylinderModel_1per.soln ../CylinderModel.fwd"
#
#
echo "#### START FORWARD LOCAL TEST AT $(date)" | tee $test_dir/level1.out
#
#
echo "#### COMMAND LINE: $cmd "  | tee -a $test_dir/level1.out
#
#
eval $cmd | tee -a $test_dir/level1.out
#
fi

#
cd ..
#
#
exit 0
#
# END OF SCRIPT

