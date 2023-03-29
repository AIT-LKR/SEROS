##########################################################################
## Makefile.
##
## The present Makefile is a pure configuration file, in which
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the Palabos user's guide for instructions on manual
## compilation.
##########################################################################

# USE: multiple arguments are separated by spaces.
#   For example: projectFiles = file1.cpp file2.cpp
#                optimFlags   = -O -finline-functions

# Leading directory of the Palabos source code
palabosRoot  = ../palabos
# Name of source files in current directory to compile and link with Palabos
SimulationFile?= seros3D.cpp
# Additional files to be copied to cluster
HeaderFiles?= ../tools/tools.h* ../tools/plb_tools.h* ../tools/plb_physicalFlowParam.h
easyppmFiles= ../libraries/easyppm/*
# Set optimization flags on/off
optimize?= true
# Set debug mode and debug flags on/off
debug?= true
# Set profiling flags on/off
profile      = false
# Set MPI-parallel mode on/off (parallelism in cluster-like environment)
MPIparallel  = true
# Set SMP-parallel mode on/off (shared-memory parallelism)
SMPparallel  = false
# Decide whether to include calls to the POSIX API. On non-POSIX systems,
#   including Windows, this flag must be false, unless a POSIX environment is
#   emulated (such as with Cygwin).
usePOSIX     = true

# Path to external source files (other than Palabos)
srcPaths =
# Path to external libraries (other than Palabos)
libraryPaths =
# Path to inlude directories (other than Palabos)
includePaths =
# Dynamic and static libraries (other than Palabos)
libraries    = -lgsl -lgslcblas -lm
# install Gnu Scientific Library "libgsl-dev" via apt

# Compiler to use without MPI parallelism
serialCXX    = g++
# Compiler to use with MPI parallelism
parallelCXX  = mpicxx
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
compileFlags = -Wall -Wnon-virtual-dtor -Wno-deprecated-declarations -std=gnu++14
# General linker flags (don't put library includes into this flag)
linkFlags    =
# Compiler flags to use when optimization mode is on
optimFlags   = -O3
# Compiler flags to use when debug mode is on
debugFlags   = -g
# Compiler flags to use when profile mode is on
profileFlags = -pg


##########################################################################
# The code below is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(palabosRoot)/scons/scons.py -j 6 -f $(palabosRoot)/SConstruct

SConsArgs = palabosRoot=$(palabosRoot) \
            projectFiles="$(SimulationFile)" \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            SMPparallel=$(SMPparallel) \
            usePOSIX=$(usePOSIX) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
            profileFlags="$(profileFlags)" \
            srcPaths="$(srcPaths)" \
            libraryPaths="$(libraryPaths)" \
            includePaths="$(includePaths)" \
            libraries="$(libraries)"

compile:
	@echo "Compiling projectFiles..."
	python2 $(SCons) $(SConsArgs)

clean:
	python2 $(SCons) -c $(SConsArgs)
	/bin/rm -vf `find $(palabosRoot) -name '*~'`

##########################################################################
# The code below is about copying to the cluster and compiling remotely
##########################################################################

TRANSFER=scp
REXEC=ssh
SSHFLAGS=-C # Compress data
REMOTE=prak2@192.168.25.
IP?=111
DIR=~/work/praktikanten/BlacherD/bin

remote : $(SimulationFile)
	@echo "Safe-Copying projectFiles..."
	$(TRANSFER) $(SSHFLAGS) $(HeaderFiles) $(REMOTE)$(IP):$(DIR)/tools
	$(TRANSFER) $(SSHFLAGS) $(easyppmFiles) $(REMOTE)$(IP):$(DIR)/easyppm
	$(TRANSFER) $(SSHFLAGS) $(SimulationFile) $(REMOTE)$(IP):$(DIR)/pipe3D
	@echo ""; echo "Compiling projectFiles..."
	$(REXEC) $(SSHFLAGS) $(REMOTE)$(IP) "cd $(DIR) && make CompileFiles=pipe3D/$(SimulationFile)"
