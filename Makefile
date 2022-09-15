
BASEDIR=/tier2/home/groups/rhi/STAR/software
FASTJETDIR=${BASEDIR}/fastjet-install
_FJCONTRIB=${FASTJETDIR}/include/fastjet/contrib:${FJCONTRIB}
STARPICOPATH=${BASEDIR}/eventStructuredAu

AN_setter=${AN_COMMON}/AN-common-config
io_setter=${IO_LIB}/iolib-config

ccflg=`${FASTJET3}/fastjet-config --cxxflags` `root-config --cflags` `${io_setter} -I` `${AN_setter} -I`  -I${ROOUNFOLD}/src \
	  -I$(STARPICOPATH) -I$(_FJCONTRIB)/RecursiveTools

LIB_FASTJET=`${FASTJET3}/fastjet-config --cxxflags --libs`
LIB_ROOT=`root-config --cflags --glibs`
LIB_TRI= ${LIB_ROOT} ${LIB_FASTJET} `${io_setter} -L` `${AN_setter} -L` -L${ROOUNFOLD} -lRooUnfold -L$(FASTJETDIR)/lib -L${STARPICOPATH} -L$(_FJCONTRIB) -lTStarJetPico

# compilation option
CC=g++
CFLAGS=-std=c++11 -O3 -Wno-deprecated
CFLAGS_CHECK=-std=c++11 -O0 -Wno-deprecated -g

bin/runQA: obj/runQA.o obj/pAuFunctions.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI} 

bin/track_check: obj/track_check.o obj/pAuFunctions.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI} 

bin/trig_check: obj/trig_check.o obj/pAuFunctions.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI} 

bin/bbc_check_v2: obj/bbc_check_v2.o obj/pAuFunctions.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI} 

check:
	echo beans
	echo -L${STARPICOPATH}
	echo ${LIB_TRI}
	



clean:
	rm obj/* bin/*

obj/runQA.o: src/runQA.cxx src/pAu_params.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/track_check.o: src/track_check.cxx src/pAu_params.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/trig_check.o: src/trig_check.cxx src/pAu_params.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/bbc_check_v2.o: src/bbc_check_v2.cxx src/pAu_params.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/pAuFunctions.o: src/pAuFunctions.cxx src/pAuFunctions.hh
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@
