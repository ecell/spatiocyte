abs_top_builddir=@abs_top_builddir@

libdir=@libdir@

SUBDIRS = scripting RandomLib-1.6/src . tests

ECELL_DIRNAME=@ECELL_DIRNAME@
ECELL_INCLUDE_DIR=$(includedir)/$(ECELL_DIRNAME)

GSL_LIBADD = -lgsl -lgslcblas
DMTOOL_INCLUDE = -I@DMTOOL_INCLUDE_DIR@
INCLTDL = @INCLTDL@
LIBLTDL = @LIBLTDL@

RANDOMLIB = RandomLib-1.6
INCLUDES = -I$(top_srcdir)  -I$(top_srcdir)/.. -I$(srcdir)/$(RANDOMLIB)/include $(DMTOOL_INCLUDE) $(INCLTDL)
lib_LTLIBRARIES = libecs.la

AM_CXXFLAGS = $(CXXFLAGS)

libecs_la_LDFLAGS= -version-info 2:0:0 -no-undefined
libecs_la_LIBADD=  $(LIBLTDL) $(GSL_LIBADD) scripting/libecs-scripting.la RandomLib-1.6/src/libRandom.la
libecs_la_SOURCES=\
libecs.cpp\
Util.cpp\
FullID.cpp\
PropertySlot.cpp\
PropertySlotProxy.cpp\
PropertyInterface.cpp\
EcsObject.cpp\
Variable.cpp\
Stepper.cpp\
DifferentialStepper.cpp\
AdaptiveDifferentialStepper.cpp\
DiscreteTimeStepper.cpp\
DiscreteEventStepper.cpp\
PassiveStepper.cpp\
SystemStepper.cpp\
Interpolant.cpp\
Model.cpp\
EntityType.cpp\
System.cpp\
Entity.cpp\
Process.cpp\
ContinuousProcess.cpp\
DataPoint.cpp\
DataPointVector.cpp\
VVector.cpp\
osif.cpp\
Logger.cpp\
LoggerAdapter.cpp\
PhysicalLogger.cpp\
VariableReference.cpp\
LoggerBroker.cpp\
Exceptions.cpp\
StepperEvent.cpp\
SpatiocyteVector.cpp\
SpatiocyteSpecies.cpp\
SpatiocyteThread.cpp\
SpatiocyteStepper.cpp\
SpatiocyteProcess.cpp\
VisualizationLogProcess.cpp\
ReactionProcess.cpp\
DiffusionInfluencedReactionProcess.cpp\
MoleculePopulateProcess.cpp\
DiffusionProcess.cpp

RANDOMLIB_INC=$(RANDOMLIB)/include/RandomLib
randomdir=$(ECELL_INCLUDE_DIR)/RandomLib
random_HEADERS=\
$(RANDOMLIB_INC)/Random.hpp\
$(RANDOMLIB_INC)/RandomCanonical.hpp\
$(RANDOMLIB_INC)/RandomPower2.hpp\
$(RANDOMLIB_INC)/RandomEngine.hpp\
$(RANDOMLIB_INC)/RandomAlgorithm.hpp\
$(RANDOMLIB_INC)/RandomMixer.hpp\
$(RANDOMLIB_INC)/RandomSeed.hpp\
$(RANDOMLIB_INC)/RandomType.hpp

libecsdir=$(ECELL_INCLUDE_DIR)/libecs
libecs_HEADERS=\
AssocVector.h\
libecs.hpp\
Defs.hpp\
convertTo.hpp\
Polymorph.hpp\
FullID.hpp\
PropertySlot.hpp\
PropertySlotProxy.hpp\
PropertyInterface.hpp\
EcsObject.hpp\
Variable.hpp\
Stepper.hpp\
ProcessEvent.hpp\
StepperEvent.hpp\
EventScheduler.hpp\
DifferentialStepper.hpp\
AdaptiveDifferentialStepper.hpp\
DiscreteTimeStepper.hpp\
DiscreteEventStepper.hpp\
PassiveStepper.hpp\
SystemStepper.hpp\
Interpolant.hpp\
Model.hpp\
EntityType.hpp\
System.hpp\
Entity.hpp\
Process.hpp\
ContinuousProcess.hpp\
DataPoint.hpp\
DataPointVector.hpp\
PhysicalLogger.hpp\
Logger.hpp\
LoggerAdapter.hpp\
LoggerBroker.hpp\
Exceptions.hpp\
Util.hpp\
VariableReference.hpp\
DynamicPriorityQueue.hpp\
MethodProxy.hpp\
EcsObjectMaker.hpp\
PropertySlotProxyLoggerAdapter.hpp\
PropertyAttributes.hpp\
Handle.hpp\
RealMath.hpp\
SpatiocyteVector.hpp\
SpatiocyteSpecies.hpp\
SpatiocyteThread.hpp\
SpatiocyteStepper.hpp\
SpatiocyteProcess.hpp\
VisualizationLogProcess.hpp\
ReactionProcess.hpp\
DiffusionInfluencedReactionProcess.hpp\
MoleculePopulateProcess.hpp\
DiffusionProcess.hpp


noinst_HEADERS=\
osif.h\
VVector.h

ECELL_STAGING_HOME=$(abs_top_builddir)/tmp
OUTPUT=$(ECELL_STAGING_HOME)$(libdir)/libecs$(LTDL_SHLIB_EXT)
all-local: $(OUTPUT)

$(OUTPUT): libecs.la
	test -d $(ECELL_STAGING_HOME)$(libdir) || mkdir -p -- $(ECELL_STAGING_HOME)$(libdir)
	$(LIBTOOL) --mode=install cp $< $(ECELL_STAGING_HOME)$(libdir)
