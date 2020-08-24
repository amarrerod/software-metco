
if ExistsVMPI
MPICC=${mpi}
parallel = metcoPar
else
if Exists2VMPI
MPICC=${mpi}
parallel = metcoPar
else
MPICC=mpiCC
parallel = metcoPar
endif
endif

if condPlugins
temp=${plugin_path}
else
temp=${HOME}/metco
endif
includedir = ${temp}


bin_PROGRAMS= metcoSeq $(parallel)
nodist_include_HEADERS = *.so
dist-hook:
	mkdir -p $(distdir)/src/skeleton
	mkdir -p $(distdir)/src/skeleton/main
	mv $(distdir)/*.cpp $(distdir)/src/skeleton/main
	mv $(distdir)/configure* $(distdir)/src/skeleton/main
	mv $(distdir)/Makefile* $(distdir)/src/skeleton/main
	mv $(distdir)/aclocal.m4 $(distdir)/src/skeleton/main
	make -f makefile_copyplugins DIST_DIR=$(distdir)
	cp ../*.cpp $(distdir)/src/skeleton
	cp ../*.h $(distdir)/src/skeleton
	cp ../*.c $(distdir)/src/skeleton
	cp ../*.y $(distdir)/src/skeleton
	cp ../*.l $(distdir)/src/skeleton
	cp ./makefile_* $(distdir)/src/skeleton/main
	mkdir -p $(distdir)/src/skeleton/utils/
	cp ../utils/*.h ../utils/*.cpp $(distdir)/src/skeleton/utils/
	cp -R ../../../src_external $(distdir)/


bindir = ${prefix}

#ACLOCAL_AMFLAGS = -I configureFiles

metcoSeq_SOURCES = metcoSeq.cpp
metcoPar_SOURCES = metcoPar.cpp
##main_SOURCES =

OBJDIR = objdir
SODIR = sodir
LIBRARYDYR = library
FLEX = flex
BISON = bison
CCFLAGS = -c -O2 -g -std=c++14 
IFLAGS = -I. -I.. -I../utils -I../../../src_external/hv-1.1-src
DFLAGS =
LFLAGS = -fPIC

all: dirs sequential ${parallel} problems algorithms conversions scores migrationSelectors exchangeSelectors mutations crossovers migrationTopologies initPopIslandLoaders localSearches outputPrinters multiObjectivizations decompositions

*.so:
../../problems/*:

dirs: objdir sodir library

library: 
	@echo Generando $@
	mkdir -p $(LIBRARYDYR)

objdir: 
	@echo Generando $@
	mkdir -p $(OBJDIR)

sodir:
	@echo Generando $@
	mkdir -p $(SODIR)

sequential: metcoSeq
parallel: metcoPar
problems:
	make -s -f makefile_problems CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
algorithms:
	make -s -f makefile_algorithms CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
conversions:
	make -s -f makefile_conversions CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
scores:
	make -s -f makefile_scores CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
migrationSelectors:
	make -s -f makefile_migration_selectors CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
exchangeSelectors:
	make -s -f makefile_exchange_selectors CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
mutations:
	make -s -f makefile_mutations CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
localSearches:
	make -s -f makefile_localSearches CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
crossovers:
	make -s -f makefile_crossovers CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
migrationTopologies:
	make -s -f makefile_migration_topologies CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
initPopIslandLoaders:
	make -s -f makefile_init_pop_island_loaders CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
outputPrinters:
	make -s -f makefile_output_printers CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
multiObjectivizations:
	make -s -f makefile_multiobjectivizations CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"
decompositions:
		make -s -f makefile_decompositions CXX="$(CXX)" CCFLAGS="$(CCFLAGS)"

.PHONY: all clean


######################## sequential #####################################

metcoSeq: metcoSeq.cpp $(OBJDIR)/Individual.o $(OBJDIR)/MOFrontBinaryInteger.o $(OBJDIR)/MOFrontVector_noMPI.o $(OBJDIR)/MOFront.o $(OBJDIR)/hv.o $(OBJDIR)/avl.o $(OBJDIR)/EA.o $(OBJDIR)/utils.o $(OBJDIR)/ScoreAlgorithm.o $(OBJDIR)/Mutation.o $(OBJDIR)/Crossover.o $(OBJDIR)/LocalSearch.o $(OBJDIR)/Plugin.o $(OBJDIR)/OutputPrinter.o $(OBJDIR)/MultiObjectivization.o $(OBJDIR)/Decomposition.o
	$(CXX) -DNO_CONF_PARSING $(IFLAGS) -O2 -g -std=c++14 metcoSeq.cpp $(OBJDIR)/Individual.o $(OBJDIR)/MOFrontBinaryInteger.o $(OBJDIR)/MOFrontVector_noMPI.o $(OBJDIR)/MOFront.o $(OBJDIR)/hv.o $(OBJDIR)/avl.o $(OBJDIR)/EA.o $(OBJDIR)/utils.o $(OBJDIR)/ScoreAlgorithm.o $(OBJDIR)/Mutation.o $(OBJDIR)/Crossover.o $(OBJDIR)/LocalSearch.o $(OBJDIR)/Plugin.o $(OBJDIR)/OutputPrinter.o $(OBJDIR)/MultiObjectivization.o $(OBJDIR)/Decomposition.o -o $@ -ldl



#########################################################################

############################# parallel ##################################

metcoPar: metcoPar.cpp $(OBJDIR)/MOFrontBinaryInteger.o $(OBJDIR)/MOFrontVector.o $(OBJDIR)/MOFront.o $(OBJDIR)/hv.o $(OBJDIR)/avl.o $(OBJDIR)/Individual.o $(OBJDIR)/EA.o $(OBJDIR)/ExecutionIsland.o $(OBJDIR)/CoordinatorIsland.o $(OBJDIR)/Configuration.o $(OBJDIR)/utils.o parseConf.tab.c lex.Conf.c $(OBJDIR)/Conversion.o $(OBJDIR)/ScoreAlgorithm.o $(OBJDIR)/GlobalScoreAlgorithm.o $(OBJDIR)/LocalScoreAlgorithm.o $(OBJDIR)/MigrationSelector.o $(OBJDIR)/ExchangeSelector.o $(OBJDIR)/Mutation.o $(OBJDIR)/Crossover.o $(OBJDIR)/LocalSearch.o $(OBJDIR)/MigrationTopology.o $(OBJDIR)/InitPopIslandLoader.o $(OBJDIR)/Plugin.o $(OBJDIR)/MultiObjectivization.o $(OBJDIR)/Decomposition.o
	$(MPICC) $(IFLAGS) metcoPar.cpp $(OBJDIR)/Individual.o $(OBJDIR)/MOFrontBinaryInteger.o $(OBJDIR)/MOFrontVector.o $(OBJDIR)/MOFront.o $(OBJDIR)/hv.o $(OBJDIR)/avl.o $(OBJDIR)/EA.o $(OBJDIR)/ExecutionIsland.o $(OBJDIR)/utils.o $(OBJDIR)/CoordinatorIsland.o $(OBJDIR)/Configuration.o parseConf.tab.c lex.Conf.c $(OBJDIR)/Conversion.o $(OBJDIR)/ScoreAlgorithm.o $(OBJDIR)/MigrationSelector.o $(OBJDIR)/ExchangeSelector.o $(OBJDIR)/GlobalScoreAlgorithm.o $(OBJDIR)/LocalScoreAlgorithm.o $(OBJDIR)/Mutation.o $(OBJDIR)/Crossover.o $(OBJDIR)/LocalSearch.o $(OBJDIR)/MigrationTopology.o $(OBJDIR)/InitPopIslandLoader.o $(OBJDIR)/Plugin.o $(OBJDIR)/MultiObjectivization.o $(OBJDIR)/Decomposition.o -o $@ -ldl


#########################################################################

############################# skeleton ##################################

parseConf.tab.c: ../parseConf.y
	$(BISON) -p Conf -d ../parseConf.y

lex.Conf.c: ../parseConf.l
	$(FLEX) -PConf -i ../parseConf.l

$(OBJDIR)/%.o: ../%.cpp
	$(CXX) $(IFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)%.reloc.o: ../%.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)/utils.o: ../utils/utils.cpp
	$(CXX) $(IFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)/CoordinatorIsland.o: ../CoordinatorIsland.cpp
	$(MPICC) $(IFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)/ExecutionIsland.o: ../ExecutionIsland.cpp
	$(MPICC) $(IFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)/MOFront.o: ../MOFront.cpp
	$(MPICC) $(IFLAGS) $(CCFLAGS) -DVARIANT=4 $< -o $@

$(OBJDIR)/MOFrontVector_noMPI.o: ../MOFrontVector.cpp
	$(CXX) $(IFLAGS) $(CCFLAGS) -DNO_MPI -DVARIANT=4 $< -o $@

$(OBJDIR)/MOFrontVector.o: ../MOFrontVector.cpp
	$(MPICC) $(IFLAGS) $(CCFLAGS) -DVARIANT=4 $< -o $@

$(OBJDIR)/MOFrontBinaryInteger.o: ../MOFrontBinaryInteger.cpp
	$(MPICC) $(IFLAGS) $(CCFLAGS) -DVARIANT=4 $< -o $@


$(OBJDIR)/%.reloc.o: ../%.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)/MOFrontVector_noMPI.reloc.o: ../MOFrontVector.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) -DNO_MPI -DVARIANT=4 $< -o $@

$(OBJDIR)/MOFront.reloc.o: ../MOFront.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) -DVARIANT=4 $< -o $@

$(OBJDIR)/MOFrontBinaryInteger.reloc.o: ../MOFrontBinaryInteger.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) -DVARIANT=4 $< -o $@

$(OBJDIR)/hv.reloc.o: ../../../src_external/hv-1.1-src/hv.c
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) -DVARIANT=4  $< -o $@

$(OBJDIR)/avl.reloc.o: ../../../src_external/hv-1.1-src/avl.c
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) $< -o $@

$(OBJDIR)/utils.reloc.o: ../utils/utils.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CCFLAGS) $< -o $@

#########################################################################

## Hipervolumen

$(OBJDIR)/hv.o: ../../../src_external/hv-1.1-src/hv.c
	$(CXX) $(IFLAGS) $(CCFLAGS) -DVARIANT=4 $< -o $@

$(OBJDIR)/avl.o: ../../../src_external/hv-1.1-src/avl.c
	$(CXX) $(IFLAGS) $(CCFLAGS) $< -o $@

#########################################################################
clean:
	@rm -rf $(OBJDIR)
	@rm -rf $(SODIR)
	@rm -rf $(LIBRARYDYR)
	@rm -f *.stackdump *.tmp *~ out metcoSeq metcoPar parseConf.tab.h parseConf.tab.c lex.Conf.c
	make -s -f makefile_algorithms clean
	make -s -f makefile_conversions clean
	make -s -f makefile_crossovers clean
	make -s -f makefile_exchange_selectors clean
	make -s -f makefile_migration_selectors clean
	make -s -f makefile_mutations clean
	make -s -f makefile_localSearches clean
	make -s -f makefile_problems clean
	make -s -f makefile_scores clean
	make -s -f makefile_output_printers clean
	make -s -f makefile_multiobjectivizations clean
	make -s -f makefile_decompositions clean

# Dependencias creadas usando makedepend
metcoPar.o: ../CoordinatorIsland.h ../Configuration.h ../MigrationSelector.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h ../ExchangeSelector.h ../ScoreAlgorithm.h ../Conversion.h ../EA.h ../utils/utils.h ../MOFrontVector.h ../MOFrontBinaryInteger.h ../MOFront.h ../LocalScoreAlgorithm.h ../ExecutionIsland.h
metcoSeq.o: ../EA.h ../Plugin.h ../utils/utils.h ../Individual.h ../Definitions.h ../Mutation.h ../Crossover.h ../MOFrontVector.h ../MOFrontBinaryInteger.h ../MOFront.h ../LocalScoreAlgorithm.h ../ScoreAlgorithm.h
CoordinatorIsland.o: ../CoordinatorIsland.h ../Configuration.h ../MigrationSelector.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h ../ExchangeSelector.h ../ScoreAlgorithm.h ../Conversion.h ../EA.h ../utils/utils.h ../MOFrontVector.h ../MOFrontBinaryInteger.h ../MOFront.h ../LocalScoreAlgorithm.h ../ExecutionIsland.h ../GlobalScoreAlgorithm.h
Configuration.o: ../Configuration.h ../MigrationSelector.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h ../ExchangeSelector.h ../ScoreAlgorithm.h ../Conversion.h ../EA.h ../utils/utils.h ../MOFrontBinaryInteger.h ../MOFrontVector.h ../MOFront.h ../LocalScoreAlgorithm.h
ExecutionIsland.o: ../Configuration.h ../MigrationSelector.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h ../ExchangeSelector.h ../ScoreAlgorithm.h ../Conversion.h ../EA.h ../utils/utils.h ../MOFrontBinaryInteger.h ../MOFrontVector.h ../MOFront.h ../LocalScoreAlgorithm.h ../ExecutionIsland.h ../CoordinatorIsland.h
EA.o: ../EA.h ../Plugin.h ../utils/utils.h ../Individual.h ../Definitions.h ../Mutation.h ../Crossover.h ../MOFrontBinaryInteger.h ../MOFrontVector.h ../MOFront.h ../LocalScoreAlgorithm.h ../ScoreAlgorithm.h ../Configuration.h ../MigrationSelector.h ../ExchangeSelector.h ../Conversion.h ../LocalSearch.h
MOFrontVector.o: ../MOFrontVector.h ../MOFront.h ../Definitions.h ../Individual.h ../Plugin.h ../Mutation.h ../Crossover.h ../../../src_external/hv-1.1-src/hv.h ../utils/utils.h
MOFrontBinaryInteger.o: ../MOFrontBinaryInteger.h ../MOFront.h ../Definitions.h ../Individual.h ../Plugin.h ../Mutation.h ../Crossover.h ../../../src_external/hv-1.1-src/hv.h ../utils/utils.h

MOFront.o: ../MOFront.h ../Definitions.h ../Individual.h ../Plugin.h ../Mutation.h ../Crossover.h ../../../src_external/hv-1.1-src/hv.h ../utils/utils.h
Individual.o: ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h
Conversion.o: ../Conversion.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h
ScoreAlgorithm.o: ../ScoreAlgorithm.h ../Plugin.h
GlobalScoreAlgorithm.o: ../GlobalScoreAlgorithm.h ../ScoreAlgorithm.h ../Plugin.h ../MOFront.h ../Definitions.h ../Individual.h ../Mutation.h ../Crossover.h
LocalScoreAlgorithm.o: ../LocalScoreAlgorithm.h ../ScoreAlgorithm.h ../Plugin.h ../MOFront.h ../Definitions.h ../Individual.h ../Mutation.h ../Crossover.h
MigrationSelector.o: ../MigrationSelector.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h
ExchangeSelector.o: ../ExchangeSelector.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h
Mutation.o: ../Mutation.h ../Individual.h ../Plugin.h ../Definitions.h ../Crossover.h
Crossover.o: ../Crossover.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h
parseConf.tab.o: ../Conversion.h ../Individual.h ../Plugin.h ../Definitions.h ../Mutation.h ../Crossover.h ../ScoreAlgorithm.h ../MigrationSelector.h ../ExchangeSelector.h ../Configuration.h ../EA.h ../utils/utils.h ../MOFront.h ../MOFrontVector.h ../MOFrontBinaryInteger.h ../LocalScoreAlgorithm.h
hv.o: ../../../src_external/hv-1.1-src/hv.h ../../../src_external/hv-1.1-src/avl.h
avl.o: ../../../src_external/hv-1.1-src/avl.h
Plugin.o: ../Plugin.h
