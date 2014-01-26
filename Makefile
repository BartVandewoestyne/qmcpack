include make.inc

all: qmclib tests timings

.PHONY: qmclib tests timings clean realclean

qmclib:
	cd src; ${MAKE}

tests:
	cd tests; ${MAKE}

timings:
	cd timings; ${MAKE}

clean:
	cd src; ${MAKE} clean
	cd tests; ${MAKE} clean
	cd timings; ${MAKE} clean

realclean:
	cd src; ${MAKE} realclean
	cd tests; ${MAKE} realclean
	cd timings; ${MAKE} realclean

help:
	@echo 'To build only the library:'
	@echo
	@echo '		make qmclib'
	@echo
	@echo 'To build the tests:'
	@echo
	@echo '		make tests'
	@echo
	@echo 'To build the timings:'
	@echo
	@echo '		make timings'
	@echo
	@echo 'To clean up:'
	@echo
	@echo '		make clean'
	@echo
	@echo 'To really clean up:'
	@echo
	@echo '		make realclean'
