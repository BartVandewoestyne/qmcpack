include ../make.inc

# There are currently two kinds of tests:
#
#   - test_*: tests that simply test different routines from the library.
#   - pub_*: tests that check if we can reproduce results from certain articles.
test_programs:=$(basename $(wildcard test*.f95))
article_programs:=$(basename $(wildcard [pub|write]*.f95))


all: $(test_programs) $(article_programs)

testrun: $(test_programs)
	$(foreach test, $(test_programs),./$(test);)

articletest: $(article_programs)
	$(foreach test, $(article_programs),./$(test);)

$(basename %.f95): %.f95 ../src/lib$(QMCPACKLIB).a
	$(FC) -I../src -L../src -o $@ $(FFLAGS) $< -l$(QMCPACKLIB)


.PHONY: clean
clean:
	-rm -f *.o *.mod *.dat depend.inc

realclean: clean
	-rm -f $(test_programs) $(article_programs) *.txt *.tex
