
#all: repr.pdf hecke.pdf spectrum.pdf spectrum2.pdf spectrum3.pdf qec2019.pdf

#out: repr
#	open repr.pdf 

spectrum2_Quantum.pdf: spectrum2_Quantum.tex refs3.bib 
	pdflatex spectrum2_Quantum.tex
	bibtex spectrum2_Quantum
	pdflatex spectrum2_Quantum.tex
	pdflatex spectrum2_Quantum.tex


spectrum2_IOP.pdf: spectrum2_IOP.tex refs3.bib 
	pdflatex spectrum2_IOP.tex
	bibtex spectrum2_IOP
	pdflatex spectrum2_IOP.tex
	pdflatex spectrum2_IOP.tex


qec2019.pdf: qec2019.tex refs3.bib 
	pdflatex qec2019.tex
	bibtex qec2019
	pdflatex qec2019.tex
	pdflatex qec2019.tex


spectrum3.pdf: spectrum3.tex refs3.bib pic-gap.pdf pic-gap-stabs.pdf
	pdflatex spectrum3.tex
	bibtex spectrum3
	pdflatex spectrum3.tex
	pdflatex spectrum3.tex


spectrum2.pdf: spectrum2.tex refs3.bib pic-gap.pdf pic-gap-stabs.pdf
	rm -f spectrum2.bbl
	rm -f spectrum2.aux
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

