include ../make.inc

time_programs:=$(basename $(wildcard time*.f95))


all: $(time_programs)

timing: $(time_programs)
	$(foreach timing, $(timing_programs),./$(timing);)

$(basename %.f95): %.f95 ../src/lib$(QMCPACKLIB).a
	$(FC) -I../src -L../src -o $@ $(FFLAGS) $< -l$(QMCPACKLIB)


.PHONY: clean realclean
clean:
	-rm -f *.o *.mod *.out depend.inc

realclean: clean
	-rm -f $(time_programs) *.dat
