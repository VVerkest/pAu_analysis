os = $(shell uname -s)

INCFLAGS      = -I$(shell root-config --incdir) $(shell fastjet-config --cxxflags) -I$(STARPICOPATH) -I/opt/local/include # -I$(FJCONTRIB)/RecursiveTools

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
LIBPATH       = -L$(FASTJETDIR)/lib -L$(STARPICOPATH) $(shell root-config --libs) #-L$(FJCONTRIB)
LIBS          =  $(ROOTLIBS) $(FJLIBS) -lfastjet -lfastjettools -lTStarJetPico # -I$(FJCONTRIB)

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
all : $(BDIR)/pAuHTjetUE


$(ODIR)/pAuFunctions.o : $(SDIR)/pAuFunctions.cxx $(SDIR)/pAuFunctions.hh
$(ODIR)/pAuHTjetUE.o : $(SDIR)/pAuHTjetUE.cxx $(SDIR)/pAuParameters.hh

#data analysis
$(BDIR)/pAuHTjetUE :	$(ODIR)/pAuHTjetUE.o	$(ODIR)/pAuFunctions.o



###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*
