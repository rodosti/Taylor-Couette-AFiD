#Makefile for parallel compiling
#=====================================================
# Note: use -C=all to catch array out of bounds errors
#============================================================================ 
# Compiler options
#============================================================================
# IFORT
#FC = h5pfc -r8 -ip -ipo -O3 -fpp -xAVX 


# GNU
FC = h5pfc -O3 -cpp -fdefault-real-8 -fdefault-double-8 -fopenmp 
#FC = h5pfc -O0 -cpp -fdefault-real-8 -fdefault-double-8 -fcheck=bounds

# CRAY

# GENERAL
FC += -DDOUBLE_PREC
#FC += -DDEBUG
#FC += -fopenmp 

#=======================================================================
# Library
#======================================================================

#LINKS = -L/share/apps/fftw-3.3.4/lib -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lpthread -lmkl_core -lm
LINKS = -llapack -lblas -lfftw3

#============================================================================ 
# make PROGRAM   
#============================================================================

PROGRAM = boutnp

OBJECTS = CalcMaxCFL.o AuxiliaryRoutines.o MpiAuxRoutines.o HdfRoutines.o \
          CalcLocalDivergence.o CheckDivergence.o ReadInputFile.o \
          ExplicitTermsQT.o ExplicitTermsQR.o ExplicitTermsQZ.o \
          ImplicitAndUpdateQT.o ImplicitAndUpdateQR.o ImplicitAndUpdateQZ.o \
          ReadFlowField.o CreateInitialConditions.o DumpThCuts.o InitTimeMarchScheme.o \
          CreateGrid.o OpenLogs.o StatReadReduceWrite.o \
          main.o param.o CorrectPressure.o SolveImpEqnUpdate_QT.o SolveImpEqnUpdate_QR.o \
          SolveImpEqnUpdate_QZ.o StatRoutines.o \
          HdfReadContinua.o DeallocateVariables.o WriteGridInfo.o \
          TimeMarcher.o CorrectVelocity.o CalcVelocityStats.o InitPressureSolver.o \
          SolvePressureCorrection.o CalcDissipation.o interp.o LocateLargeDivergence.o  CalcVorticity.o \
 	    decomp_2d.o decomp_2d_fft.o WriteFlowField.o CalcTorque.o QuitRoutine.o DumpFlowField.o \
 	    DumpRCuts.o  SlipBCs.o InnerCylVelocity.o InitVariables.o DumpStreamwiseAvg.o
#          alloc.o decomp_2d.o fft_fftw3.o halo.o halo_common.o

MODULES = param.o decomp_2d.o decomp_2d_fft.o
#============================================================================ 
# Linking    
#============================================================================

$(PROGRAM) : $(MODULES) $(OBJECTS)
	$(FC) $(OP_LINK) $(OBJECTS) -o $@ $(LINKS)  

#============================================================================
#  Dependencies
#============================================================================

param.o: param.F90
	$(FC) -c $(OP_COMP) param.F90

decomp_2d.o: decomp_2d.f90
	$(FC) -c $(OP_COMP) decomp_2d.f90

decomp_2d_fft.o: decomp_2d_fft.f90
	$(FC) -c $(OP_COMP) decomp_2d_fft.f90


%.o:   %.F $(MODULES)
	$(FC) -c $(OP_COMP) $< 

%.o:   %.f90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 

%.o:   %.F90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 

%.o: imbdry_codes/%.F $(MODULES)
	$(FC) -c $(OP_COMP) $< 
        

#============================================================================
#  Clean up
#============================================================================

clean :
	rm *.o  
	rm *.mod

veryclean :
	rm *.o *.mod *.out *.h5 stats/*.h5 dati/* boutnp
