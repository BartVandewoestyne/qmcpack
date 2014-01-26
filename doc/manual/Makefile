# Adopt these to your own needs...
LATEX = latex
BIBTEX = bibtex
PDFLATEX = ps2pdf
DVIPS = dvips

# The targets
TARGET = qmcpack_manual
SOURCES = $(wildcard *.tex)

all: dvi ps pdf

dvi: ${TARGET}.dvi

ps: ${TARGET}.ps

pdf: ${TARGET}.pdf

# Device Independent
${TARGET}.dvi: ${SOURCES}
	echo "Running latex..."
	${LATEX} ${TARGET}.tex
	echo "Running bibtex..."
	${BIBTEX} ${TARGET}
	echo "Rerunning latex...."
	${LATEX} ${TARGET}.tex
	echo "Rerunning latex...."
	${LATEX} ${TARGET}.tex

# Postscript
${TARGET}.ps: ${TARGET}.dvi
	${DVIPS} -Ppdf -o ${TARGET}.ps ${TARGET}.dvi

# Portable Document Format
${TARGET}.pdf: ${TARGET}.ps
	${PDFLATEX} ${TARGET}.ps


.PHONY: clean realclean
clean:
	rm -f *~ *.aux *.log *.bbl *.toc *.blg

realclean: clean
	rm -f *.dvi *.ps *.pdf
