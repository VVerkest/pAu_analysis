os = $(shell uname -s)

INCFLAGS      = -I$(shell root-config --incdir) $(shell fastjet-config --cxxflags) -I$(STARPICOPATH) -I$(FJCONTRIB)/RecursiveTools -I/opt/local/include

ifeq ($(os),Linux)
CXXFLAGS      = -std=c++17
else
CXXFLAGS      = -O -std=c++11 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)
FJLIBS        = $(shell fastjet-config --plugins=yes --libs)
LIBPATH       = -L$(FASTJETDIR)/lib -L$(STARPICOPATH) $(shell root-config --libs) -L$(FJCONTRIB)
LIBS          =  $(ROOTLIBS) $(FJLIBS) -I$(FJCONTRIB) -lfastjet -lfastjettools -lTStarJetPico

# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################


###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o
	@echo
	@echo LINKING
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBPATH) $(LIBS)

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all : $(BDIR)/pAu_QA $(BDIR)/pAu_analysis_MB $(BDIR)/pAu_analysis_HT $(BDIR)/pAu_HT_jets $(BDIR)/pAu_HT_dijets $(BDIR)/pAu_HT_jetTree $(BDIR)/pAu_HT_jetPlot $(BDIR)/find_bad_towers $(BDIR)/find_bad_trigger_towers $(BDIR)/pAuHTjetUE $(BDIR)/pAuHTdijetUE $(BDIR)/EAdistribution $(BDIR)/pAu_HT_jet_trigger

$(ODIR)/pAuQAFunctions.o : $(SDIR)/pAuQAFunctions.cxx $(SDIR)/pAuQAFunctions.hh
$(ODIR)/pAu_QA.o : $(SDIR)/pAu_QA.cxx $(SDIR)/pAuQAFunctions.hh
$(ODIR)/pAuFunctions.o : $(SDIR)/pAuFunctions.cxx $(SDIR)/pAuFunctions.hh
$(ODIR)/pAu_analysis_MB.o : $(SDIR)/pAu_analysis_MB.cxx
$(ODIR)/pAu_analysis_HT.o : $(SDIR)/pAu_analysis_HT.cxx
$(ODIR)/pAu_HT_jets.o : $(SDIR)/pAu_HT_jets.cxx $(SDIR)/pAu_HT_jetParameters.hh $(SDIR)/bemc_helper.h
$(ODIR)/pAu_HT_jet_trigger.o : $(SDIR)/pAu_HT_jet_trigger.cxx $(SDIR)/pAu_HT_jetParameters.hh
$(ODIR)/pAuHTjetUE.o : $(SDIR)/pAuHTjetUE.cxx $(SDIR)/pAu_HT_jetParameters.hh
$(ODIR)/pAuHTdijetUE.o : $(SDIR)/pAuHTdijetUE.cxx $(SDIR)/pAu_HT_jetParameters.hh
$(ODIR)/EAdistribution.o : $(SDIR)/EAdistribution.cxx $(SDIR)/pAu_HT_jetParameters.hh
$(ODIR)/pAu_HT_dijets.o : $(SDIR)/pAu_HT_dijets.cxx
$(ODIR)/pAu_HT_jetTree.o : $(SDIR)/pAu_HT_jetTree.cxx
$(ODIR)/pAu_HT_jetPlot.o : $(SDIR)/pAu_HT_jetPlot.cxx
$(ODIR)/find_bad_towers.o : $(SDIR)/find_bad_towers.cxx
$(ODIR)/find_bad_trigger_towers.o : $(SDIR)/find_bad_trigger_towers.cxx

#data analysis
$(BDIR)/pAu_QA :	$(ODIR)/pAu_QA.o	$(ODIR)/pAuQAFunctions.o
$(BDIR)/pAu_analysis_MB :	$(ODIR)/pAu_analysis_MB.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAu_analysis_HT :	$(ODIR)/pAu_analysis_HT.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAu_HT_jets :	$(ODIR)/pAu_HT_jets.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAu_HT_jet_trigger :	$(ODIR)/pAu_HT_jet_trigger.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAuHTjetUE :	$(ODIR)/pAuHTjetUE.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAuHTdijetUE :	$(ODIR)/pAuHTdijetUE.o	$(ODIR)/pAuFunctions.o
$(BDIR)/EAdistribution :	$(ODIR)/EAdistribution.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAu_HT_dijets :	$(ODIR)/pAu_HT_dijets.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAu_HT_jetTree :	$(ODIR)/pAu_HT_jetTree.o	$(ODIR)/pAuFunctions.o
$(BDIR)/pAu_HT_jetPlot :	$(ODIR)/pAu_HT_jetPlot.o
$(BDIR)/find_bad_towers :	$(ODIR)/find_bad_towers.o
$(BDIR)/find_bad_trigger_towers :	$(ODIR)/find_bad_trigger_towers.o


###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*
