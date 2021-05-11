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
all : $(BDIR)/sim $(BDIR)/ptResponse_diff $(BDIR)/ptResponse_diff_OLD  $(BDIR)/ptResponse $(BDIR)/ptDifferentialUE $(BDIR)/JEScorrection

$(ODIR)/funcs.o : $(SDIR)/funcs.cxx $(SDIR)/funcs.hh $(SDIR)/params.hh
$(ODIR)/UEfuncs.o : $(SDIR)/UEfuncs.cxx $(SDIR)/UEfuncs.hh $(SDIR)/params.hh
$(ODIR)/sim.o : $(SDIR)/sim.cxx $(SDIR)/params.hh
$(ODIR)/ptResponse_diff_OLD.o : $(SDIR)/ptResponse_diff_OLD.cxx $(SDIR)/params.hh
$(ODIR)/ptResponse_diff.o : $(SDIR)/ptResponse_diff.cxx $(SDIR)/params.hh
$(ODIR)/ptResponse.o : $(SDIR)/ptResponse.cxx $(SDIR)/params.hh
$(ODIR)/ptDifferentialUE.o : $(SDIR)/ptDifferentialUE.cxx $(SDIR)/params.hh
$(ODIR)/JEScorrection.o : $(SDIR)/JEScorrection.cxx $(SDIR)/params.hh

#data analysis
$(BDIR)/sim :	$(ODIR)/sim.o	$(ODIR)/funcs.o
$(BDIR)/ptResponse_diff_OLD :       $(ODIR)/ptResponse_diff_OLD.o       $(ODIR)/funcs.o
$(BDIR)/ptResponse_diff :	$(ODIR)/ptResponse_diff.o	$(ODIR)/funcs.o
$(BDIR)/ptResponse :	$(ODIR)/ptResponse.o	$(ODIR)/funcs.o
$(BDIR)/ptDifferentialUE :	$(ODIR)/ptDifferentialUE.o	$(ODIR)/UEfuncs.o
$(BDIR)/JEScorrection :	$(ODIR)/JEScorrection.o	$(ODIR)/funcs.o


###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*
