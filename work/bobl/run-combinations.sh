cd /home/ellen.porter/cnc-framework/examples/Combinations

[ -r Makefile ] || ucnc_t -p icnc/tcp

make run WORKLOAD_ARGS="1000 500"
