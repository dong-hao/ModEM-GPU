# Makefile suited for building the Mod3DMT program
# Generated using: ./fmkmf.pl [OPTIONS] Mod3DMT.f90 > Makefile
# with command line options
# -p .:MPI:INV:LAPACK:SENS:UTILS:FIELDS/FiniteDiff3D:3D_MT/FWD_SP:3D_MT/SP_Topology:3D_MT:3D_MT/DICT:3D_MT/modelParam:3D_MT/FWD:3D_MT/FWD/Mod2d:3D_MT/ioMod
# -f90 ifort (compiler)
# -opt -O3 -w  -xSSE4.2 -std03 -i8 -I/opt/intel/lib (compiler optimisation)
# -lp /opt/intel/mkl/lib (linking options: path to libraries)
# -l -lmkl_lapack95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core  -lpthread -qopenmp -parallel (linking options)
# -o ./objs/3D_MT/IFortReleaseMPI_SP (output directory for object files)

#  Uncomment these lines to make program for Solaris OS (legacy)
# F90 = f90
# FFLAGS = -dalign -g -C -w  -L/usr/local/lib
# LIBS = -xlic_lib=sunperf
#  Uncomment these lines to make program with g95
# include Makefile.local
# OBJDIR = ./objs/3D_MT/G95Debug
# F90 = g95
# FFLAGS = -O2
# FFLAGS = -g -ftrace=frame -fbounds-check
# MPIFLAGS = -cpp # for serial code
# MODULE = -fmod=$(OBJDIR)
# LIBS = -lblas -llapack
#  Uncomment these lines to make program with Intel compiler
# include Makefile.local
# OBJDIR = ./objs/3D_MT/IFortDebug
# F90 = ifort
# FFLAGS = -O3 -parallel -openmp #-heap-arrays
# FFLAGS = -debug all -check bounds -traceback -heap-arrays
# MPIFLAGS = -cpp # for serial code
# MODULE = -module $(OBJDIR)
# LIBS = -lblas -llapack
#  Uncomment these lines to make program with PGI compiler
# include Makefile.local
# OBJDIR = ./objs/3D_MT/PGIDebug
# F90 = pgf95  # mpif90
# FFLAGS = -O3
# FFLAGS = -g -Mprof=lines -Mbounds
# MPIFLAGS = -Mpreprocess # for serial code
# MPIFLAGS = -Bstatic  -Mipa=fast  -Mextend  -Kieee -Mpreprocess -DMPI
# MODULE = -module $(OBJDIR)
# LIBS = -llapack -lblas
# LIBS = -L/usr/lib64 -llapack -lblas -lpgftnrtl -Mprof=lines

# ------------------Macro-Defs---------------------
include Makefile.local
OBJDIR = ./objs/3D_MT/GFortReleaseMPI_SP
F90 =  mpif90
FFLAGS = -O3 -ffree-line-length-none -mavx
MPIFLAGS = -x f95-cpp-input -DMPI
MODULE = -J $(OBJDIR)
LIBS_PATH = -L/usr/lib -L/usr/local/opt/lapack/lib
LIBS = -llapack -lblas

# -------------------End-macro-Defs---------------------------
OBJ = $(OBJDIR)/math_constants.o $(OBJDIR)/fields_orientation.o $(OBJDIR)/utilities.o $(OBJDIR)/file_units.o $(OBJDIR)/DataSpace.o $(OBJDIR)/GridDef.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_spherical.o $(OBJDIR)/elements.o $(OBJDIR)/GridCalc.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/Declaration_MPI.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/transmitters.o $(OBJDIR)/receivers.o $(OBJDIR)/SensMatrix.o $(OBJDIR)/EMfieldInterp.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/dataTypes.o $(OBJDIR)/DataFunc.o $(OBJDIR)/DataSens.o $(OBJDIR)/SolverSens.o $(OBJDIR)/vecTranslate.o $(OBJDIR)/spOpTools.o $(OBJDIR)/spOpTopology_SG.o $(OBJDIR)/MetricElements_CSG.o $(OBJDIR)/nestedEM.o $(OBJDIR)/WSfwd2Dpar.o $(OBJDIR)/WSutils.o $(OBJDIR)/WSfwd1Dmod.o $(OBJDIR)/WSfwd2Dmod.o $(OBJDIR)/FwdTEmod.o $(OBJDIR)/boundary_ws.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/solver.o $(OBJDIR)/EMsolve3D.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/SensComp.o $(OBJDIR)/UserCtrl.o $(OBJDIR)/Sub_MPI.o $(OBJDIR)/Main_MPI.o $(OBJDIR)/SymmetryTest.o $(OBJDIR)/ioAscii.o $(OBJDIR)/DataIO.o $(OBJDIR)/Main.o $(OBJDIR)/INVcore.o $(OBJDIR)/NLCG.o $(OBJDIR)/DCG.o $(OBJDIR)/LBFGS.o $(OBJDIR)/Mod3DMT_SP.o 


all: Mod3DMT 

# Here is the link step 
Mod3DMT: $(OBJDIR) $(OBJ) 
	 $(F90) -o $(OUTDIR)/Mod3DMT_SP $(OBJ) $(LIBS_PATH) $(LIBS)

# Here are the compile steps 

$(OBJDIR): 
	mkdir -p $(OBJDIR)

$(OBJDIR)/math_constants.o:UTILS/math_constants.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/math_constants.f90 -o $(OBJDIR)/math_constants.o

$(OBJDIR)/fields_orientation.o:UTILS/fields_orientation.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/fields_orientation.f90 -o $(OBJDIR)/fields_orientation.o

$(OBJDIR)/utilities.o:UTILS/utilities.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/utilities.f90 -o $(OBJDIR)/utilities.o

$(OBJDIR)/file_units.o:UTILS/file_units.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/file_units.f90 -o $(OBJDIR)/file_units.o

$(OBJDIR)/DataSpace.o:SENS/DataSpace.f90 $(OBJDIR)/utilities.o $(OBJDIR)/math_constants.o $(OBJDIR)/fields_orientation.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SENS/DataSpace.f90 -o $(OBJDIR)/DataSpace.o

$(OBJDIR)/GridDef.o:3D_MT/GridDef.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/GridDef.f90 -o $(OBJDIR)/GridDef.o

$(OBJDIR)/sg_vector.o:FIELDS/FiniteDiff3D/sg_vector.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/GridDef.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) FIELDS/FiniteDiff3D/sg_vector.f90 -o $(OBJDIR)/sg_vector.o

$(OBJDIR)/sg_scalar.o:FIELDS/FiniteDiff3D/sg_scalar.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/GridDef.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) FIELDS/FiniteDiff3D/sg_scalar.f90 -o $(OBJDIR)/sg_scalar.o

$(OBJDIR)/sg_spherical.o:FIELDS/FiniteDiff3D/sg_spherical.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) FIELDS/FiniteDiff3D/sg_spherical.f90 -o $(OBJDIR)/sg_spherical.o

$(OBJDIR)/elements.o:UTILS/elements.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/elements.f90 -o $(OBJDIR)/elements.o

$(OBJDIR)/GridCalc.o:3D_MT/GridCalc.f90 $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_spherical.o $(OBJDIR)/elements.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/GridCalc.f90 -o $(OBJDIR)/GridCalc.o

$(OBJDIR)/sg_sparse_vector.o:FIELDS/FiniteDiff3D/sg_sparse_vector.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/GridDef.o $(OBJDIR)/sg_vector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) FIELDS/FiniteDiff3D/sg_sparse_vector.f90 -o $(OBJDIR)/sg_sparse_vector.o

$(OBJDIR)/Declaration_MPI.o:MPI/Declaration_MPI.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) MPI/Declaration_MPI.f90 -o $(OBJDIR)/Declaration_MPI.o

$(OBJDIR)/ModelSpace.o:3D_MT/modelParam/ModelSpace.f90 $(OBJDIR)/GridCalc.o $(OBJDIR)/file_units.o $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/Declaration_MPI.o 3D_MT/modelParam/modelCov/RecursiveAR.hd 3D_MT/modelParam/ModelMap.inc 3D_MT/modelParam/modelCov/RecursiveAR.inc 3D_MT/modelParam/modelParamIO/Binary.inc 3D_MT/modelParam/modelParamIO/Mackie.inc 3D_MT/modelParam/modelParamIO/WS.inc 3D_MT/modelParam/ModelParam_MPI.inc
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/modelParam/ModelSpace.f90 -o $(OBJDIR)/ModelSpace.o

$(OBJDIR)/transmitters.o:3D_MT/DICT/transmitters.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/DICT/transmitters.f90 -o $(OBJDIR)/transmitters.o

$(OBJDIR)/receivers.o:3D_MT/DICT/receivers.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/DICT/receivers.f90 -o $(OBJDIR)/receivers.o

$(OBJDIR)/SensMatrix.o:SENS/SensMatrix.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/file_units.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSpace.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/transmitters.o $(OBJDIR)/receivers.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SENS/SensMatrix.f90 -o $(OBJDIR)/SensMatrix.o

$(OBJDIR)/EMfieldInterp.o:3D_MT/EMfieldInterp.f90 $(OBJDIR)/utilities.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/GridCalc.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/EMfieldInterp.f90 -o $(OBJDIR)/EMfieldInterp.o

$(OBJDIR)/sg_boundary.o:FIELDS/FiniteDiff3D/sg_boundary.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/GridDef.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_sparse_vector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) FIELDS/FiniteDiff3D/sg_boundary.f90 -o $(OBJDIR)/sg_boundary.o

$(OBJDIR)/SolnSpace.o:3D_MT/SolnSpace.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/transmitters.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/SolnSpace.f90 -o $(OBJDIR)/SolnSpace.o

$(OBJDIR)/dataTypes.o:3D_MT/DICT/dataTypes.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/DICT/dataTypes.f90 -o $(OBJDIR)/dataTypes.o

$(OBJDIR)/DataFunc.o:3D_MT/DataFunc.f90 $(OBJDIR)/EMfieldInterp.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/receivers.o $(OBJDIR)/transmitters.o $(OBJDIR)/dataTypes.o $(OBJDIR)/fields_orientation.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/DataFunc.f90 -o $(OBJDIR)/DataFunc.o

$(OBJDIR)/DataSens.o:SENS/DataSens.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/fields_orientation.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSpace.o $(OBJDIR)/DataFunc.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SENS/DataSens.f90 -o $(OBJDIR)/DataSens.o

$(OBJDIR)/SolverSens.o:3D_MT/SolverSens.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/transmitters.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/SolverSens.f90 -o $(OBJDIR)/SolverSens.o

$(OBJDIR)/vecTranslate.o:3D_MT/SP_Topology/vecTranslate.f90 $(OBJDIR)/utilities.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/SP_Topology/vecTranslate.f90 -o $(OBJDIR)/vecTranslate.o

$(OBJDIR)/spOpTools.o:3D_MT/SP_Topology/spOpTools.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/SP_Topology/spOpTools.f90 -o $(OBJDIR)/spOpTools.o

$(OBJDIR)/spOpTopology_SG.o:3D_MT/SP_Topology/spOpTopology_SG.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/SP_Topology/spOpTopology_SG.f90 -o $(OBJDIR)/spOpTopology_SG.o

$(OBJDIR)/MetricElements_CSG.o:3D_MT/SP_Topology/MetricElements_CSG.f90 $(OBJDIR)/GridCalc.o $(OBJDIR)/vecTranslate.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/SP_Topology/MetricElements_CSG.f90 -o $(OBJDIR)/MetricElements_CSG.o

$(OBJDIR)/nestedEM.o:3D_MT/FWD/nestedEM.f90 $(OBJDIR)/EMfieldInterp.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/sg_sparse_vector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/nestedEM.f90 -o $(OBJDIR)/nestedEM.o

$(OBJDIR)/WSfwd2Dpar.o:3D_MT/FWD/Mod2d/WSfwd2Dpar.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/Mod2d/WSfwd2Dpar.f90 -o $(OBJDIR)/WSfwd2Dpar.o

$(OBJDIR)/WSutils.o:3D_MT/FWD/Mod2d/WSutils.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/Mod2d/WSutils.f90 -o $(OBJDIR)/WSutils.o

$(OBJDIR)/WSfwd1Dmod.o:3D_MT/FWD/Mod2d/WSfwd1Dmod.f90 $(OBJDIR)/WSfwd2Dpar.o $(OBJDIR)/WSutils.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/Mod2d/WSfwd1Dmod.f90 -o $(OBJDIR)/WSfwd1Dmod.o

$(OBJDIR)/WSfwd2Dmod.o:3D_MT/FWD/Mod2d/WSfwd2Dmod.f90 $(OBJDIR)/WSfwd1Dmod.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/Mod2d/WSfwd2Dmod.f90 -o $(OBJDIR)/WSfwd2Dmod.o

$(OBJDIR)/FwdTEmod.o:3D_MT/FWD/Mod2d/FwdTEmod.f90 $(OBJDIR)/WSfwd2Dmod.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/Mod2d/FwdTEmod.f90 -o $(OBJDIR)/FwdTEmod.o

$(OBJDIR)/boundary_ws.o:3D_MT/FWD/boundary_ws.f90 $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/FwdTEmod.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD/boundary_ws.f90 -o $(OBJDIR)/boundary_ws.o

$(OBJDIR)/modelOperator3D.o:3D_MT/FWD_SP/modelOperator3D.f90 $(OBJDIR)/vecTranslate.o $(OBJDIR)/spOpTools.o $(OBJDIR)/spOpTopology_SG.o $(OBJDIR)/MetricElements_CSG.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/nestedEM.o $(OBJDIR)/boundary_ws.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD_SP/modelOperator3D.f90 -o $(OBJDIR)/modelOperator3D.o

$(OBJDIR)/solver.o:3D_MT/FWD_SP/solver.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/spOpTools.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/modelOperator3D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD_SP/solver.f90 -o $(OBJDIR)/solver.o

$(OBJDIR)/EMsolve3D.o:3D_MT/FWD_SP/EMsolve3D.f90 $(OBJDIR)/sg_boundary.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/vecTranslate.o $(OBJDIR)/solver.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/sg_vector.o $(OBJDIR)/vecTranslate.o $(OBJDIR)/solver.o $(OBJDIR)/spOpTools.o $(OBJDIR)/modelOperator3D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/FWD_SP/EMsolve3D.f90 -o $(OBJDIR)/EMsolve3D.o

$(OBJDIR)/ForwardSolver.o:3D_MT/ForwardSolver.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/DataFunc.o $(OBJDIR)/DataSpace.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/EMsolve3D.o $(OBJDIR)/transmitters.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/ForwardSolver.f90 -o $(OBJDIR)/ForwardSolver.o

$(OBJDIR)/SensComp.o:SENS/SensComp.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/fields_orientation.o $(OBJDIR)/utilities.o $(OBJDIR)/SensMatrix.o $(OBJDIR)/DataSens.o $(OBJDIR)/SolverSens.o $(OBJDIR)/ForwardSolver.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SENS/SensComp.f90 -o $(OBJDIR)/SensComp.o

$(OBJDIR)/UserCtrl.o:./UserCtrl.f90 $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) ./UserCtrl.f90 -o $(OBJDIR)/UserCtrl.o

$(OBJDIR)/Sub_MPI.o:3D_MT/Sub_MPI.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/UserCtrl.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/Declaration_MPI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/Sub_MPI.f90 -o $(OBJDIR)/Sub_MPI.o

$(OBJDIR)/Main_MPI.o:MPI/Main_MPI.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/file_units.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSens.o $(OBJDIR)/SolverSens.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/SensComp.o $(OBJDIR)/Declaration_MPI.o $(OBJDIR)/Sub_MPI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) MPI/Main_MPI.f90 -o $(OBJDIR)/Main_MPI.o

$(OBJDIR)/SymmetryTest.o:SENS/SymmetryTest.f90 $(OBJDIR)/SensComp.o $(OBJDIR)/Main_MPI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SENS/SymmetryTest.f90 -o $(OBJDIR)/SymmetryTest.o

$(OBJDIR)/ioAscii.o:3D_MT/ioMod/ioAscii.f90 $(OBJDIR)/GridDef.o $(OBJDIR)/math_constants.o $(OBJDIR)/EMsolve3D.o $(OBJDIR)/DataSpace.o $(OBJDIR)/DataFunc.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/transmitters.o $(OBJDIR)/receivers.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/ioMod/ioAscii.f90 -o $(OBJDIR)/ioAscii.o

$(OBJDIR)/DataIO.o:3D_MT/DataIO.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/file_units.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSpace.o $(OBJDIR)/GridCalc.o $(OBJDIR)/transmitters.o $(OBJDIR)/receivers.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/DataIO.f90 -o $(OBJDIR)/DataIO.o

$(OBJDIR)/Main.o:3D_MT/Main.f90 $(OBJDIR)/ModelSpace.o $(OBJDIR)/DataSpace.o $(OBJDIR)/DataFunc.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/SensMatrix.o $(OBJDIR)/UserCtrl.o $(OBJDIR)/ioAscii.o $(OBJDIR)/DataIO.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) 3D_MT/Main.f90 -o $(OBJDIR)/Main.o

$(OBJDIR)/INVcore.o:INV/INVcore.f90 $(OBJDIR)/SensComp.o $(OBJDIR)/DataIO.o $(OBJDIR)/Main_MPI.o $(OBJDIR)/Sub_MPI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) INV/INVcore.f90 -o $(OBJDIR)/INVcore.o

$(OBJDIR)/NLCG.o:INV/NLCG.f90 $(OBJDIR)/INVcore.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) INV/NLCG.f90 -o $(OBJDIR)/NLCG.o

$(OBJDIR)/DCG.o:INV/DCG.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/SensComp.o $(OBJDIR)/Main.o $(OBJDIR)/Main_MPI.o $(OBJDIR)/Sub_MPI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) INV/DCG.f90 -o $(OBJDIR)/DCG.o

$(OBJDIR)/LBFGS.o:INV/LBFGS.f90 $(OBJDIR)/INVcore.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) INV/LBFGS.f90 -o $(OBJDIR)/LBFGS.o

$(OBJDIR)/Mod3DMT_SP.o:Mod3DMT.f90 $(OBJDIR)/SensComp.o $(OBJDIR)/SymmetryTest.o $(OBJDIR)/Main.o $(OBJDIR)/NLCG.o $(OBJDIR)/DCG.o $(OBJDIR)/LBFGS.o $(OBJDIR)/Main_MPI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Mod3DMT.f90 -o $(OBJDIR)/Mod3DMT_SP.o

# Type " make clean " to get rid of all object and module files 
clean: 
	cd $(OBJDIR); \
	rm -f *~ *.o *.obj *.mod *.d *.s00 *.dbg *.stackdump \
	`find . -mindepth 1 -name "*~"` 

cleanall: clean 
	rm -f $(OUTDIR)/Mod3DMT_SP 

src: clean 
	tar cvfz $(ARCHIVE).tgz * 

  
