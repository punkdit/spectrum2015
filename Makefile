
all: repr.pdf hecke.pdf spectrum.pdf spectrum2.pdf

#out: repr
#	open repr.pdf 

spectrum2.pdf: spectrum2.tex refs3.bib
	pdflatex spectrum2.tex
	bibtex spectrum2
	pdflatex spectrum2.tex
	pdflatex spectrum2.tex


spectrum.pdf: spectrum.tex refs3.bib
	pdflatex spectrum.tex
	bibtex spectrum
	pdflatex spectrum.tex
	pdflatex spectrum.tex


repr.pdf: repr.tex refs.bib pic-gcolor-1.pdf
	pdflatex repr.tex
	bibtex repr
	pdflatex repr.tex
	pdflatex repr.tex

repr-abstract.pdf: repr-abstract.tex refs.bib 
	pdflatex repr-abstract.tex
	bibtex repr-abstract
	pdflatex repr-abstract.tex
	pdflatex repr-abstract.tex


hecke.pdf: hecke.tex 
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
	gcc -g -fPIC -shared -I/usr/include/python2.7 -o pnauty.so -O3 pnauty.c nauty/nauty.a -lefence -lpython2.7


cbracket.so: cbracket.c
	gcc -g -fPIC -shared -I/usr/include/python2.7 -o cbracket.so -O3 cbracket.c -lpython2.7


pic:
	./gauge.py size=1   transparency=0.4 R=2.0 filename=pic-gcolor-1.pdf
	./gauge.py size=1.5 transparency=0.5 R=2.5 filename=pic-gcolor-15.pdf
	./gauge.py size=2   transparency=0.4 R=3.0 filename=pic-gcolor-2.pdf # broken ??

