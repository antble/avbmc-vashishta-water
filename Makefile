#SHELL = /bin/sh
CUR_DIR := $(shell pwd)
# ====================================
#        Compiler Options
# ====================================
FC := mpif90
#FC := /opt/openmpi/bin/mpif90
#FC := mpifort
#FC := gfortran
CC := mpicc
OPTIMIZE_FLAGS := -O3
#OPTIMIZE_FLAGS += -xHost
#OPTIMIZE_FLAGS += -ipo
#OPTIMIZE_FLAGS += -no-prec-div
#OPTIMIZE_FLAGS += -prof-gen -prof-dir=$(CUR_DIR)/profiling
#OPTIMIZE_FLAGS += -prof-use -prof-dir=$(CUR_DIR)/profiling
# DETAILEDDEBUG:= -fbacktrace -fcheck=all -g -ffree-line-length-0 -Og
# DETAILEDDEBUG:= -check all -traceback -g -fpe
# DEBUGFLAGS:= -fbacktrace -fcheck=all -g
# DEBUGFLAGS += -heap-arrays 1024
# DEBUGFLAGS += $(DETAILEDDEBUG)
# DEBUGFLAGS += -check all -traceback -g
# DEBUGFLAGS += -pg 
# DEBUGFLAGS += -ffpe-trap=invalid
# DEBUGFLAGS += -Wunused-parameter 
# DEBUGFLAGS := -fimplicit-none  -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace
COMPFLAGS := $(DEBUGFLAGS) $(OPTIMIZE_FLAGS)


# ====================================
#        Directory List
# ====================================

SRC := $(CUR_DIR)/src
MODS := $(CUR_DIR)
OBJ := $(CUR_DIR)/objects
TRIAL := $(CUR_DIR)/Trials
ESUB := $(CUR_DIR)/src/EnergyFunctions
LJ_Q := $(CUR_DIR)/src/EnergyFunctions/LJ_Q
VASHISHTA := $(CUR_DIR)/src/EnergyFunctions/Vashishta
PEDONE := $(CUR_DIR)/src/EnergyFunctions/Pedone
TERSOFF := $(CUR_DIR)/src/EnergyFunctions/Tersoff
CBMC := $(CUR_DIR)/src/CBMC_Functions
SWAP := $(CUR_DIR)/src/SwapFunctions
ANALYSIS_SUB := $(CUR_DIR)/src/AnalysisFunctions
BOUNDARY := $(CUR_DIR)/src/BoundaryFunctions

#RUN_DIR := $(TRIAL)/Trial1_TransRot_Test
#RUN_DIR := $(TRIAL)/Trial2_Water_Test



# Define file extensions
.SUFFIXES:
.SUFFIXES: .f .f90 .o .mod 
# ====================================
#        Compiler specific commands
# ====================================

#MODFLAGS := -I $(MODS) -J $(MODS)
# ====================================
#        Source Files
# ====================================

MOD_SRC := $(SRC)/VariablePrecision.f90\
			$(SRC)/Common.f90 \
			$(SRC)/Units.f90 \
			$(SRC)/ForceFieldFunctions.f90\
			$(SRC)/DistanceStorage.f90
SRC_ENERGY := $(LJ_Q)/Bending_Functions.f90 \
            $(LJ_Q)/BondStretch_Functions.f90 \
            $(LJ_Q)/LJ_Electro_Functions.f90 \
            $(LJ_Q)/LJ_Electro_Functions_DistStore.f90 \
            $(LJ_Q)/Pressure_LJ_Electro_Functions.f90 \
            $(LJ_Q)/Intra_LJ_Electro_Functions.f90\
            $(LJ_Q)/Torsional_Functions.f90 \
            $(LJ_Q)/Improper_Functions.f90 \
            $(LJ_Q)/Rosen_Boltz_Functions.f90\
            $(LJ_Q)/EnergyInterfaceFunctions_LJ.f90\
            $(LJ_Q)/EnergyInterfaceFunctions_LJ_DistStore.f90\
			$(VASHISHTA)/Vashishta_Functions.f90\
			$(VASHISHTA)/Rosen_Vashishta_Functions.f90\
			$(VASHISHTA)/EnergyInterfaceFunctions_Vashishta.f90\
            $(PEDONE)/Pedone_Functions_Experimental.f90\
            $(PEDONE)/Rosen_Pedone_Functions.f90\
            $(PEDONE)/EnergyInterfaceFunctions_Pedone_Experimental.f90\
            $(TERSOFF)/Tersoff_Functions.f90\
            $(TERSOFF)/Rosen_Tersoff_Functions.f90\
            $(TERSOFF)/EnergyInterface_Tersoff.f90\
            $(ESUB)/EnergyPointers.f90
SRC_CRIT:=  $(SRC)/ClusterCriteria_Energy.f90\
            $(SRC)/ClusterCriteria_Distance.f90\
            $(SRC)/ClusterCriteria_Distance_PairStore.f90
#SRC_BOUNDARY := $(BOUNDARY)/ClusterCriteria_Energy_New.f90
#SRC_BOUNDARY := $(SRC)/ClusterCriteria_Energy.f90
SRC_BIAS := $(SRC)/UmbrellaSampling_Version2.f90\
            $(SRC)/Umbrella_Types.f90\
            $(SRC)/WHAM_Version2.f90
SRC_INPUT := $(SRC)/ScriptInput.f90\
			$(SRC)/ScriptForceField.f90\
			$(SRC)/Input_Ultility.f90 \
			# $(SRC)/ReadInput.f90 
SRC_MAIN := $(SRC)/BasicMovement.f90\
            $(SRC)/MCMove_Module.f90\
            $(SRC)/DebugFunctions.f90\
            $(SRC)/Main.f90\
			$(SRC)/RandomNew.f90\
			$(SRC)/RandomTools.f90\
			$(SRC)/OutputFunctions.f90\
			$(SRC)/UmbrellaSampling.f90\
			$(SRC)/AngleIntegration.f90\
			$(SRC)/CoordinateFunctions.f90\
			$(SRC)/SelfAdaptiveDistribution.f90
SRC_MAIN2:=  $(SRC)/ETableFunctions.f90
SRC_CBMC := $(CBMC)/CBMC.f90\
            $(CBMC)/CBMC_Initialize.f90\
            $(CBMC)/CBMC_ConfigGen.f90\
            $(CBMC)/CBMC_Utility.f90\
            $(CBMC)/CBMC_PartialRegrowth.f90\
            $(CBMC)/CBMC_Rosen_AVBMC_ConfigGen.f90
SRC_SWAP := $(SWAP)/Exchange.f90\
            $(SWAP)/SwapBoundaries.f90\
            $(SWAP)/AVBMC_EBias_Rosen.f90
SRC_ANALYSIS := $(ANALYSIS_SUB)/MiscelaniousVariables.f90\
            $(ANALYSIS_SUB)/SimplePairDistance.f90\
            $(ANALYSIS_SUB)/RadialDistributionFunction.f90\
            $(ANALYSIS_SUB)/Q3Functions.f90\
            $(ANALYSIS_SUB)/Q4Functions.f90\
            $(ANALYSIS_SUB)/Q6Functions.f90\
            $(ANALYSIS_SUB)/RadialDensity.f90\
            $(ANALYSIS_SUB)/AnalysisMain.f90
# SRC_SWAP := $(SWAP)/AVBMC_EBias.f90
# SRC_SWAP := $(SWAP)/AVBMC_Uniform.f90

SRC_COMPLETE := $(SRC_MAIN) $(SRC_SWAP) $(SRC_ANALYSIS) $(SRC_BIAS) $(SRC_INPUT) $(SRC_CRIT) $(SRC_MAIN2) $(SRC_ENERGY)

# ====================================
#        Object Files
# ====================================
OBJ_MAIN:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_MAIN))
OBJ_MAIN2:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_MAIN2))

OBJ_TEMP:=$(patsubst $(PEDONE)/%.f90,$(OBJ)/%.o,$(SRC_ENERGY))
OBJ_TEMP:=$(patsubst $(TERSOFF)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))
OBJ_TEMP:=$(patsubst $(LJ_Q)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))
OBJ_TEMP:=$(patsubst $(VASHISHTA)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))
OBJ_ENERGY:=$(patsubst $(ESUB)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))

OBJ_MOD:=$(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(MOD_SRC))
OBJ_CBMC:=$(patsubst $(CBMC)/%.f90,$(OBJ)/%.o,$(SRC_CBMC))
OBJ_INPUT:=$(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SRC_INPUT))
OBJ_BIAS:=$(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SRC_BIAS))
OBJ_SWAP:=$(patsubst $(SWAP)/%.f90,$(OBJ)/%.o,$(SRC_SWAP))
OBJ_CRIT:=$(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SRC_CRIT))
OBJ_ANALYSIS:=$(patsubst $(ANALYSIS_SUB)/%.f90,$(OBJ)/%.o,$(SRC_ANALYSIS))
OBJ_BOUNDARY:=$(patsubst $(BOUNDARY)/%.f90,$(OBJ)/%.o,$(SRC_BOUNDARY))

OBJ_COMPLETE:= $(OBJ_CRIT) $(OBJ_BOUNDARY) $(OBJ_ENERGY) $(OBJ_BIAS) $(OBJ_INPUT) $(OBJ_ANALYSIS) $(OBJ_MAIN2) $(OBJ_CBMC) $(OBJ_SWAP) $(OBJ_MAIN) 
# ====================================
#        Compile Commands
# ====================================


.f.o :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
.f90.o :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

.f.mod :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

.f90.mod :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(ESUB)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(TERSOFF)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(PEDONE)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(LJ_Q)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(VASHISHTA)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
        
$(OBJ)/%.o: $(CBMC)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
            
$(OBJ)/%.o: $(SRC_ANALYSIS)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(SWAP)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(ANALYSIS_SUB)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(SRC)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(BOUNDARY)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<


# ====================================
#        Compile Commands
# ====================================
default: startUP generalNucleation finale
debug: startUP_debug generalNucleation_debug finale
quick: startUP generalNucleation finale
neat: startUP generalNucleation removeObjects finale
clean: removeObjects removeExec finale    
    
createMods: $(MOD_SRC) 
		@echo =============================================
		@echo            Creating Module Files
		@echo =============================================		
		@echo		
		@echo  -------- Compiling VariablePrecision.f90
		@$(FC) -c $(SRC)/VariablePrecision.f90 $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/VariablePrecision.o
		@echo  -------- Compiling Units.f90
		@$(FC) -c $(SRC)/Units.f90  $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/Units.o		
		@echo  -------- Compiling ForceFieldFunctions.f90
		@$(FC) -c $(SRC)/ForceFieldFunctions.f90 $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/ForceFieldFunctions.o				
		@echo  -------- Compiling Common.f90
		@$(FC) -c $(SRC)/Common.f90  $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/Common.o
		@echo  -------- Compiling DistanceStorage.f90
		@$(FC) -c $(SRC)/DistanceStorage.f90 $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/DistanceStorage.o		
		@echo  -------- Compiling SwapBoundaries.f90
		@$(FC) -c $(SWAP)/SwapBoundaries.f90 $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/SwapBoundaries.o	
		@echo =============================================
		@echo            Creating Object Files
		@echo =============================================	
		@echo  	


	

generalNucleation: $(OBJ_MOD)  $(OBJ_COMPLETE) 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(COMPFLAGS) $(MODFLAGS)  $^ -o $@ 	
	

generalNucleation_debug: $(OBJ_MOD)  $(OBJ_COMPLETE) 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(DETAILEDDEBUG) $(MODFLAGS)  $^ -o $@ 	
		
startUP:
		@echo ==================================================================
		@echo ---------------------- Begin ---------------------------------		
		@echo Current Directory:$(CUR_DIR)		
		@echo Compiler and Flags used:	$(FC) $(COMPFLAGS) 		
		@echo		

startUP_debug:
		@echo ==================================================================
		@echo ---------------------- Begin ---------------------------------		
		@echo Current Directory:$(CUR_DIR)		
		@echo Compiler and Flags used:	$(FC) $(DETAILEDDEBUG) 		
		@echo		

finale:
		@echo
		@echo ---------------------- Finished! ---------------------------------
		@echo ==================================================================		

removeObjects:
		@echo =============================================
		@echo            Cleaning Directory
		@echo =============================================		
		@echo		
		@rm -f ./*.o ./*.mod				
		@rm -f $(SRC)/*.o $(SRC)/*.mod		
		@rm -f $(MODS)/*.o $(MODS)/*.mod			
		@rm -f $(SRC)/*/*.o $(SRC)/*/*.mod
		@rm -f $(OBJ)/*.o		

removeExec:
		@rm -f $(CUR_DIR)/generalNucleation
		@rm -f $(CUR_DIR)/generalNucleation.exe            
		@rm -f $(CUR_DIR)/generalNucleation_debug
		@rm -f $(CUR_DIR)/generalNucleation_debug.exe            



# ====================================
#        Dependencies
# ====================================
$(OBJ)/Common.o: $(OBJ)/VariablePrecision.o $(OBJ)/Units.o

$(OBJ)/ForceFieldFunctions.o: $(OBJ)/Common.o $(OBJ)/VariablePrecision.o

$(OBJ)/Units.o: $(OBJ)/VariablePrecision.o

$(OBJ)/MCMove_Module.o: $(OBJ_ENERGY) $(OBJ_CBMC) $(OBJ_SWAP) $(OBJ)/BasicMovement.o
$(OBJ)/Exchange.o: $(OBJ)/ETableFunctions.o


$(OBJ)/CBMC.o: $(OBJ)/Common.o $(OBJ_ENERGY) $(OBJ)/UmbrellaSampling_Version2.o  $(OBJ)/CBMC_Utility.o
$(OBJ)/CBMC_ConfigGen.o: $(OBJ)/CoordinateFunctions.o

$(OBJ)/ClusterCriteria_Distance.o: $(OBJ)/Common.o $(OBJ)/ForceFieldFunctions.o $(OBJ)/ClusterCriteria_Distance_PairStore.o
$(OBJ)/ClusterCriteria_Energy.o: $(OBJ)/Common.o $(OBJ)/ForceFieldFunctions.o

$(OBJ)/CoordinateFunctions.o: $(OBJ)/Common.o $(OBJ)/RandomTools.o $(OBJ)/CBMC_Utility.o

$(OBJ)/AnalysisMain.o: $(OBJ)/Q6Functions.o $(OBJ)/RadialDistributionFunction.o $(OBJ)/SimplePairDistance.o $(OBJ)/MiscelaniousVariables.o $(OBJ)/RadialDensity.o $(OBJ)/Common.o

$(OBJ)/SimplePairDistance.o: $(OBJ)/Umbrella_Types.o

$(OBJ)/EnergyPointers.o: $(OBJ)/EnergyInterfaceFunctions_Pedone_Experimental.o $(OBJ)/EnergyInterfaceFunctions_LJ_DistStore.o

$(OBJ)/SimplePairDistance.o: $(OBJ)/DistanceStorage.o $(OBJ)/MiscelaniousVariables.o 


$(OBJ)/EnergyInterfaceFunctions_LJ.o: $(OBJ)/LJ_Electro_Functions_DistStore.o $(OBJ)/Rosen_Boltz_Functions.o $(OBJ)/Bending_Functions.o $(OBJ)/BondStretch_Functions.o $(OBJ)/Torsional_Functions.o $(OBJ)/Intra_LJ_Electro_Functions.o $(OBJ)/Improper_Functions.o $(OBJ)/Rosen_Tersoff_Functions.o $(OBJ)/EnergyInterfaceFunctions_LJ_DistStore.o
$(OBJ)/EnergyInterfaceFunctions_Pedone_Experimental.o: $(OBJ)/Pedone_Functions_Experimental.o $(OBJ)/Rosen_Pedone_Functions.o
$(OBJ)/EnergyInterfaceFunctions_Vashishta.o:  $(OBJ)/Rosen_Vashishta_Functions.o
$(OBJ)/Pedone_Functions_Experimental.o: $(OBJ)/DistanceStorage.o
$(OBJ)/Tersoff_Functions.o: $(OBJ)/DistanceStorage.o
$(OBJ)/Rosen_Tersoff_Functions.o: $(OBJ)/Tersoff_Functions.o

$(OBJ)/MiscelaniousVariables.o: $(OBJ)/Common.o

$(OBJ)/Main.o: $(OBJ)/Common.o $(OBJ)/SelfAdaptiveDistribution.o

$(OBJ)/ScriptInput.o: $(OBJ)/Common.o $(OBJ)/UmbrellaSampling_Version2.o $(OBJ)/AnalysisMain.o $(OBJ)/MCMove_Module.o $(OBJ)/ScriptForceField.o
$(OBJ)/ScriptForceField.o: $(OBJ)/Common.o $(OBJ)/EnergyPointers.o 

$(OBJ)/Umbrella_Types.o: $(OBJ)/VariablePrecision.o 
$(OBJ)/UmbrellaSampling_Version2.o: $(OBJ)/Common.o $(OBJ)/SimplePairDistance.o $(OBJ)/Umbrella_Types.o $(OBJ_ANALYSIS)
$(OBJ)/WHAM_Version2.o: $(OBJ)/SwapBoundaries.o $(OBJ)/Common.o $(OBJ)/UmbrellaSampling_Version2.o

