# This is where "setup_env.sh" has to be sourced from.
TOP = ~ellen.porter

CNC_BINDIR = $(TOP)/cnc-framework/bin

# The default, /bin/sh, doesn't support "source", "pushd", and "popd".
SHELL = /bin/bash

TCP_NUM_CLIENTS=4

SETUP_ENV = \
	source $(TOP)/icnc/bin/cncvars.sh ;\
	pushd $(TOP)/cnc-framework ;\
	source setup_env.sh ;\
	popd ;\
	export TCP_NUM_CLIENTS=$(TCP_NUM_CLIENTS)

# compiles and runs
.PHONY: run
run: Makefile
	$(SETUP_ENV) ;\
	$(MAKE) -f Makefile run

# makes Makefile.icnc cnc_support/icnc/* *.c *.h icnc-tcp-start.sh
.PHONY: non-distributed
non-distributed: immaculate
	$(SETUP_ENV) ;\
	$(CNC_BINDIR)/ucnc_t -p icnc/x86

# makes Makefile.icnc-tcp (w/link) cnc_support/icnc/* *.c *.h icnc-tcp-start.sh
.PHONY: distributed
distributed: immaculate
	$(SETUP_ENV) ;\
	$(CNC_BINDIR)/ucnc_t -p icnc/tcp

# sets environment variables and starts a shell
.PHONY: shell
shell:
	$(SETUP_ENV) ;\
	$(SHELL)

# makes Makefile.icnc-tcp (w/link) cnc_support/icnc/* *.c *.h icnc-tcp-start.sh
.PHONY: list-cnc-options
list-cnc-options:
	$(SETUP_ENV) ;\
	$(CNC_BINDIR)/ucnc_t -h

# lists environment variables set or changed by $(SETUP_ENV) scripts
.PHONY: list-env
list-env:
	export >/tmp/env_vars\$\$\_before.txt ;\
	$(SETUP_ENV) ;\
	export >/tmp/env_vars\$\$\_after.txt ;\
	diff /tmp/env_vars\$\$\_before.txt /tmp/env_vars\$\$\_after.txt ;\
	rm /tmp/env_vars\$\$\_before.txt /tmp/env_vars\$\$\_after.txt

.PHONY: clean
clean:
	rm -f *~ core

immaculate: clean
	rm -rf Makefile* *.c *.h cnc_support/ icnc-tcp-start.sh install/ build/

