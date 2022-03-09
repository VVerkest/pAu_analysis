#!/usr/bin/bash

# ROOT
module load root

# recommended paths for use:
export BASEDIR=/tier2/home/groups/rhi/STAR/software

# FASTJETDIR
export FASTJETDIR=${BASEDIR}/fastjet-install

### FastJet Contrib
export FJCONTRIB=${FASTJETDIR}/include/fastjet/contrib:${FJCONTRIB}
export LD_LIBRARY_PATH=${FJCONTRIB}:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${FJCONTRIB}:${DYLD_LIBRARY_PATH}
export PATH=${PATH}:${FJCONTRIB}

### PYTHIA8
export PYTHIA8DIR=${BASEDIR}/pythia8
export PYTHIA8DATA=${PYTHIA8DIR}/xmldoc

export STARPICOPATH=${BASEDIR}/eventStructuredAu

###### Update paths
if [[ -n ${LD_LIBRARY_PATH} ]]; then
    export LD_LIBRARY_PATH
fi

export PATH=~/physics/analysis/jetmass2/macros/:${PATH}
export PATH=${PATH}:${FASTJETDIR}/bin:${PYTHIA8DIR}/bin
export PATH=${PATH}:${FASTJETDIR}/include/fastjet/contrib/
export PATH="${PATH}:${ROOUNFOLDDIR}/bin"
export LD_LIBRARY_PATH=${FASTJETDIR}/lib:${PYTHIA8DIR}/lib:${STARPICOPATH}:${LD_LIBRARY_PATH}

# print statement if interactive shell
if [[ -n $PS1 ]] ; then

    echo "
    Setup ROOT, ktJet (including FastJet)
    =====================================
    Setting up the following environments: 
    ROOT:         $ROOTSYS
    PYTHIA8:      $PYTHIA8DIR
    FastJet:      $FASTJETDIR
    STARPICOPATH: $STARPICOPATH"
    if [[ -n $ROOUNFOLDDIR ]]; then
        echo ROOUNFOLDDIR: $ROOUNFOLDDIR
    fi
    echo ""

fi
