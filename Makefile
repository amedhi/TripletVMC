#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
#include ../root_dir.mk # must be included first
include ./make_options.mk
#-------------------------------------------------------------
# Source files
SRCS = scheduler/mpi_comm.cpp 
SRCS+= scheduler/cmdargs.cpp 
SRCS+= scheduler/inputparams.cpp 
SRCS+= scheduler/taskparams.cpp 
SRCS+= scheduler/worker.cpp 
SRCS+= scheduler/master_scheduler.cpp
SRCS+= scheduler/scheduler.cpp
#SRCS+= xml/pugixml.cpp 
#SRCS+= expression/expression.cpp 
#SRCS+= expression/tokens.cpp 
#SRCS+= expression/functions.cpp
#SRCS+= expression/objects.cpp
SRCS+= expression/complex_expression.cpp
SRCS+= utils/utils.cpp 
SRCS+= lattice/lattice.cpp
SRCS+= lattice/latticelibrary.cpp
SRCS+= lattice/graph.cpp
SRCS+= model/hamiltonian_term.cpp
SRCS+= model/model.cpp
SRCS+= model/modellibrary.cpp
SRCS+= wavefunction/varparm.cpp
SRCS+= wavefunction/blochbasis.cpp
SRCS+= wavefunction/mf_model.cpp
SRCS+= wavefunction/groundstate.cpp
SRCS+= wavefunction/fermisea.cpp
SRCS+= wavefunction/bcs_state.cpp
SRCS+= wavefunction/wavefunction.cpp
SRCS+= wavefunction/projector.cpp
SRCS+= mcdata/mcdata.cpp
SRCS+= mcdata/mc_observable.cpp
SRCS+= vmc/random.cpp
SRCS+= vmc/basisstate.cpp
SRCS+= vmc/sysconfig.cpp
SRCS+= vmc/energy.cpp
SRCS+= vmc/sccorr.cpp
SRCS+= vmc/particle.cpp
SRCS+= vmc/observables.cpp
SRCS+= vmc/stochastic_reconf.cpp
SRCS+= optimizer/optimizer.cpp
#SRCS+= vmc/measurement.cpp
SRCS+= vmc/vmc.cpp
SRCS+= vmc/simulator.cpp
SRCS+= main.cpp
VMC_SRCS = $(addprefix src/,$(SRCS))
#-------------------------------------------------------------
# Headers
HDRS=    scheduler/mpi_comm.h \
         scheduler/optionparser.h scheduler/cmdargs.h \
         scheduler/inputparams.h scheduler/worker.h scheduler/task.h \
         scheduler/scheduler.h \
         expression/complex_expression.h \
         utils/utils.h \
         lattice/constants.h lattice/lattice.h lattice/graph.h \
	 montecarlo/simulator.h \
         model/modelparams.h  model/quantum_op.h \
	 model/hamiltonian_term.h \
	 model/model.h \
	 wavefunction/varparm.h \
	 wavefunction/blochbasis.h \
	 wavefunction/mf_model.h \
	 wavefunction/groundstate.h \
	 wavefunction/fermisea.h \
	 wavefunction/bcs_state.h \
	 wavefunction/disordered_sc.h \
	 variational/wavefunction.h \
	 wavefunction/projector.h \
	 mcdata/mcdata.h  \
	 mcdata/mc_observable.h  \
	 vmc/energy.h \
	 vmc/sccorr.h \
	 vmc/particle.h \
	 vmc/observables.h \
	 vmc/random.h  vmc/basisstate.h vmc/sysconfig.h \
	 vmc/disorder.h \
	 vmc/stochastic_reconf.h \
	 vmc/vmc.h \
	 vmc/simulator.h \
	 vmcpp.h \
#         model/qn.h model/quantum_operator.h model/sitebasis.h model/modelparams.h \
#         observable/mcdata.h observable/observables.h \
#	 montecarlo/observable_operator.h montecarlo/simulator.h \
	 simulation.h
VMC_HDRS = $(addprefix src/,$(HDRS))
MUPARSER_LIB = $(PROJECT_ROOT)/src/expression/muparserx/libmuparserx.a
#-------------------------------------------------------------
# Target
TAGT=a.out
ifeq ($(MPI), HAVE_BOOST_MPI)
TAGT=v.out
endif

# Put all auto generated stuff to this build dir.
ifeq ($(BUILD_DIR), $(CURDIR))
  $(error In-source build is not allowed, choose another build directory)
endif

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cpp,$(BUILD_DIR)/%.o,$(VMC_SRCS))
# GCC/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 
# compiler flags

.PHONY: all
all: $(TAGT) #$(INCL_HDRS)

$(TAGT): $(OBJS) $(MUPARSER_LIB)
	$(VMC_CXX) -o $(TAGT) $(OBJS) $(VMC_LDFLAGS) $(VMC_LIBS)  $(MUPARSER_LIB)

%.o: %.cpp
	$(VMC_CXX) -c $(VMC_CXXFLAGS) -o $@ $<

# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	@echo "$(VMC_CXX) -c $(VMC_CXXFLAGS) -o $(@F) $(<F)"
	@$(VMC_CXX) -MMD -c $(VMC_CXXFLAGS) -o $@ $<

$(VMC_INCLDIR)/%.h: %.h 
	@mkdir -p $(@D)
	@echo "Copying $< to 'include'" 
	@cp -f $< $@

$(MUPARSER_LIB):
	@cd ./src/expression/muparserx/ && $(MAKE)

# installation
#prefix = ../install#/usr/local
#libdir = $(prefix)/lib
#includedir = $(prefix)/include/cmc++

.PHONY: install
install:	
	@echo "Already installed in $(VMC_LIBDIR) and $(VMC_INCLDIR)" 

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TAGT)"
	@rm -f $(TAGT) 

.PHONY: bclean
bclean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
