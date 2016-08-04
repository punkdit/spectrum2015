
all: repr.pdf hecke.pdf

#out: repr
#	open repr.pdf 

repr.pdf: repr.tex refs.bib
	pdflatex repr.tex
	bibtex repr
	pdflatex repr.tex
	pdflatex repr.tex


hecke.pdf: hecke.tex refs.bib
	pdflatex hecke.tex
	bibtex hecke
	pdflatex hecke.tex
	pdflatex hecke.tex


notes: notes.tex refs.bib
	pdflatex notes.tex
	bibtex notes
	pdflatex notes.tex
	pdflatex notes.tex





spectrum: spectrum.tex refs.bib
	pdflatex spectrum.tex
	bibtex spectrum
	pdflatex spectrum.tex
	pdflatex spectrum.tex
	open spectrum.pdf 

pnauty.so: pnauty.c
	gcc -g -shared -I/usr/include/python2.7 -o pnauty.so -O3 -march=i686 pnauty.c nauty/nauty.a -lefence -lpython2.7


