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
all : $(BDIR)/JEScorrection $(BDIR)/newJEScorrection $(BDIR)/noJEScorrection $(BDIR)/noJEScorrection_prelimTrackEffic

$(ODIR)/funcs.o : $(SDIR)/funcs.cxx $(SDIR)/funcs.hh $(SDIR)/params.hh
$(ODIR)/JEScorrection.o : $(SDIR)/JEScorrection.cxx $(SDIR)/params.hh
$(ODIR)/newJEScorrection.o : $(SDIR)/newJEScorrection.cxx $(SDIR)/params.hh
$(ODIR)/noJEScorrection.o : $(SDIR)/noJEScorrection.cxx $(SDIR)/params.hh
$(ODIR)/noJEScorrection_prelimTrackEffic.o : $(SDIR)/noJEScorrection_prelimTrackEffic.cxx $(SDIR)/params.hh

#data analysis
$(BDIR)/JEScorrection :	$(ODIR)/JEScorrection.o	$(ODIR)/funcs.o
$(BDIR)/newJEScorrection :	$(ODIR)/newJEScorrection.o	$(ODIR)/funcs.o
$(BDIR)/noJEScorrection :	$(ODIR)/noJEScorrection.o	$(ODIR)/funcs.o
$(BDIR)/noJEScorrection_prelimTrackEffic :	$(ODIR)/noJEScorrection_prelimTrackEffic.o	$(ODIR)/funcs.o

###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*
