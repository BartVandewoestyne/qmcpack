include ../make.inc

FOBJ:=$(patsubst %.f95,%.o,$(wildcard *.f95))

all: qmclib install

qmclib: depend $(FOBJ)
	$(ARCH) $(ARCHFLAGS) lib$(QMCPACKLIB).a $(FOBJ)
	$(RANLIB) lib$(QMCPACKLIB).a

install:
	mkdir -p $(INSTALLDIR)/qmcpack
	cp lib$(QMCPACKLIB).a $(INSTALLDIR)/qmcpack
	cp *.mod $(INSTALLDIR)/qmcpack

include depend.inc

%.o: %.f95
	$(FC) -c $(FFLAGS) $<

.PHONY: clean realclean depend

depend: depend.inc

depend.inc:
	makedepf90 *.f95 > depend.inc

clean:
	-rm -f *.o *.mod

realclean: clean
	-rm -f depend.inc *.a
