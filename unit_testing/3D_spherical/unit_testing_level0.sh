#!/bin/bash
#
# Forward modeling setup was created using the unit_testing_setup.m script in Matlab. Cartesian and
# spherical models were setup to mimic each other, however, there is no meaningful conversion between
# the two. However, the results should be directly comparable.
#
# Comparison with ModEMM done using FWD_3D_spherical.m script in Matlab.
#
# Unit testing level 0: runs the symmetry tests and the derivative test. Should run and produce self-
# consistent outputs. Does not require comparison against other outputs.

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

# IF NEEDED, UNCOMMENT TO REGENERATE INPUT FILES
# ../$EXEC -F CylinderModel+QuarterSpace+Spherical.rho CylinderModel_6per_Spherical_Template.dat CylinderModel_6per.dat CylinderModel_6per.soln CylinderModel.fwd
# ../$EXEC -F CylinderModel+QuarterSpace+Spherical.rho CylinderModel_1per_Spherical_Template.dat CylinderModel_1per.dat CylinderModel_1per.soln CylinderModel.fwd
#
# CREATE TEST OUTPUT FOLDER
mkdir -p OUTPUT
mkdir -p "OUTPUT/$test_dir"
#
# ENTER TEST OUTPUT FOLDER
cd "OUTPUT"

if [ $# -gt 2 ]; then
#
# PREPARE FOR TESTING
../$EXEC  -A  m ../CylinderModel+QuarterSpace+Spherical.rho ../delta.rho 0.02
../$EXEC  -A  d ../CylinderModel_6per.dat ../delta.dat 0.02
../$EXEC  -A  e ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_6per.dat ../CylinderModel_6per.soln ../delta.soln 0.02
../$EXEC  -b  ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_6per.dat ../efield.bc ../CylinderModel.fwd
../$EXEC  -A  b ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_6per.dat ../efield.bc ../delta.bc 0.02 ../CylinderModel.fwd
# 
#
# SET UP FILE NAMES
rFile_Model=../CylinderModel+QuarterSpace+Spherical.rho
rFile_Data=../CylinderModel_6per.dat
rFile_EMsoln=../CylinderModel_6per.soln
rFile_dEMsoln=../delta.soln
rFile_EMrhs=../delta.bc
rFile_dModel=../delta.rho
rFile_fwdCtrl=../CylinderModel.fwd
#
#
echo "#### START SYMMETRY MPI TESTS WITH $ncores CORES AT $(date) ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: mpirun -n $ncores ../$EXEC  -A  L $rFile_Model $rFile_dEMsoln $rFile_Data $test_dir/LTd.soln $test_dir/Le.dat $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
# 
mpirun -n $ncores ../$EXEC  -A  L $rFile_Model $rFile_dEMsoln $rFile_Data $test_dir/LTd.soln $test_dir/Le.dat $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality d^T L e = e^T L^T d for any EMsoln and data.
#   Optionally, outputs L e and L^T d.
#
#
echo "#### COMMAND LINE: mpirun -n $ncores ../$EXEC  -A  S $rFile_Model $rFile_EMrhs $rFile_Data $test_dir/Sb.soln $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#  
mpirun -n $ncores ../$EXEC  -A  S $rFile_Model $rFile_EMrhs $rFile_Data $test_dir/Sb.soln $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality b^T S^{-1} b = b^T (S^{-1})^T b for any EMrhs.
#   For simplicity, use one EMrhs for forward and transpose solvers.
#   Data file only needed to set up dictionaries.
#   Optionally, outputs e = S^{-1} b.
#
#
echo "#### COMMAND LINE: mpirun -n $ncores ../$EXEC  -A  P $rFile_Model $rFile_dModel $rFile_EMsoln $rFile_Data $test_dir/PTe.rho $test_dir/Pm.rhs $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#  
mpirun -n $ncores ../$EXEC  -A  P $rFile_Model $rFile_dModel $rFile_EMsoln $rFile_Data $test_dir/PTe.rho $test_dir/Pm.rhs $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality e^T P m = m^T P^T e for any EMsoln and data.
#   The data template isn't needed here except to set up the transmitters.
#   Optionally, outputs P m and P^T e.
#
#
echo "#### COMMAND LINE: mpirun -n $ncores ../$EXEC  -A  Q $rFile_Model $rFile_dModel $rFile_Data $test_dir/QTd.rho $test_dir/Qm.dat $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#  
mpirun -n $ncores ../$EXEC  -A  Q $rFile_Model $rFile_dModel $rFile_Data $test_dir/QTd.rho $test_dir/Qm.dat $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality d^T Q m = m^T Q^T d for any model and data.
#   Optionally, outputs Q m and Q^T d.
#
#
echo "#### COMMAND LINE: mpirun -n $ncores ../$EXEC  -A  O $rFile_Model $rFile_Data $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
# 
mpirun -n $ncores ../$EXEC  -A  O $rFile_Model $rFile_Data $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests all intermediate operators: grad, curl, div and grid elements.
#
#
echo "#### COMMAND LINE: mpirun -n $ncores ../$EXEC  -A  J $rFile_Model $rFile_dModel $rFile_Data $test_dir/JTd.rho $test_dir/Jm.dat $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#
mpirun -n $ncores ../$EXEC  -A  J $rFile_Model $rFile_dModel $rFile_Data $test_dir/JTd.rho $test_dir/Jm.dat $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality d^T J m = m^T J^T d for any model and data.
#   Optionally, outputs J m and J^T d.#
#
else
#
# PREPARE FOR TESTING
../$EXEC  -A  m ../CylinderModel+QuarterSpace+Spherical.rho ../delta.rho 0.02
../$EXEC  -A  d ../CylinderModel_1per.dat ../delta.dat 0.02
../$EXEC  -A  e ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_1per.dat ../CylinderModel_1per.soln ../delta.soln 0.02
../$EXEC  -b  ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_1per.dat ../efield.bc ../CylinderModel.fwd
../$EXEC  -A  b ../CylinderModel+QuarterSpace+Spherical.rho ../CylinderModel_1per.dat ../efield.bc ../delta.bc 0.02 ../CylinderModel.fwd
# 
#
# SET UP FILE NAMES
rFile_Model=../CylinderModel+QuarterSpace+Spherical.rho
rFile_Data=../CylinderModel_1per.dat
rFile_EMsoln=../CylinderModel_1per.soln
rFile_dEMsoln=../delta.soln
rFile_EMrhs=../delta.bc
rFile_dModel=../delta.rho
rFile_fwdCtrl=../CylinderModel.fwd
#
#
echo "#### START SYMMETRY LOCAL TESTS AT $(date) ####" | tee $test_dir/level0.out
#
#
echo "#### COMMAND LINE: ../$EXEC  -A  L $rFile_Model $rFile_dEMsoln $rFile_Data $test_dir/LTd.soln $test_dir/Le.dat $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
# 
../$EXEC  -A  L $rFile_Model $rFile_dEMsoln $rFile_Data $test_dir/LTd.soln $test_dir/Le.dat $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality d^T L e = e^T L^T d for any EMsoln and data.
#   Optionally, outputs L e and L^T d.
#
#
echo "#### COMMAND LINE: ../$EXEC  -A  S $rFile_Model $rFile_EMrhs $rFile_Data $test_dir/Sb.soln $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#  
../$EXEC  -A  S $rFile_Model $rFile_EMrhs $rFile_Data $test_dir/Sb.soln $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality b^T S^{-1} b = b^T (S^{-1})^T b for any EMrhs.
#   For simplicity, use one EMrhs for forward and transpose solvers.
#   Data file only needed to set up dictionaries.
#   Optionally, outputs e = S^{-1} b.
#
#
echo "#### COMMAND LINE: ../$EXEC  -A  P $rFile_Model $rFile_dModel $rFile_EMsoln $rFile_Data $test_dir/PTe.rho $test_dir/Pm.rhs $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#  
../$EXEC  -A  P $rFile_Model $rFile_dModel $rFile_EMsoln $rFile_Data $test_dir/PTe.rho $test_dir/Pm.rhs $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality e^T P m = m^T P^T e for any EMsoln and data.
#   The data template isn't needed here except to set up the transmitters.
#   Optionally, outputs P m and P^T e.
#
#
echo "#### COMMAND LINE: ../$EXEC  -A  Q $rFile_Model $rFile_dModel $rFile_Data $test_dir/QTd.rho $test_dir/Qm.dat $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#  
../$EXEC  -A  Q $rFile_Model $rFile_dModel $rFile_Data $test_dir/QTd.rho $test_dir/Qm.dat $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality d^T Q m = m^T Q^T d for any model and data.
#   Optionally, outputs Q m and Q^T d.
#
#
echo "#### COMMAND LINE: ../$EXEC  -A  O $rFile_Model $rFile_Data $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
# 
../$EXEC  -A  O $rFile_Model $rFile_Data $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests all intermediate operators: grad, curl, div and grid elements.
#
#
#
echo "#### COMMAND LINE: ../$EXEC  -A  J $rFile_Model $rFile_dModel $rFile_Data $test_dir/JTd.rho $test_dir/Jm.dat $rFile_fwdCtrl"  | tee -a $test_dir/level0.out
#
../$EXEC  -A  J $rFile_Model $rFile_dModel $rFile_Data $test_dir/JTd.rho $test_dir/Jm.dat $rFile_fwdCtrl | tee -a $test_dir/level0.out
#   Tests the equality d^T J m = m^T J^T d for any model and data.
#   Optionally, outputs J m and J^T d.
#
fi
#
#   Generate quick summary output
#
../../symmetry_tests_summary.sh $test_dir/level0.out
#
#
cd ..
#
exit 0
#
# END OF SCRIPT
