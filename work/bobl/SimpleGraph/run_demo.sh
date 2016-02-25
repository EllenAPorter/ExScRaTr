#! /bin/bash

top=~ellen.porter

source $top/icnc/bin/cncvars.sh

#echo $CNCROOT
#echo $TBBROOT

pushd $top/cnc-framework
source setup_env.sh
popd

export TCP_NUM_CLIENTS=4
echo "TCP_NUM_CLIENTS:" $TCP_NUM_CLIENTS

#echo $UCNC_ROOT

# non-distributed
# make Makefile build/ cnc_support/ install/
#ucnc_t -p icnc

# compiles and runs
#make run

# distributed
# make Makefile build/ cnc_support/ install/
ucnc_t -p icnc/tcp

# compiles and runs
make run

