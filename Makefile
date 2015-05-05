

out: repr
	open repr.pdf 

repr: repr.tex refs.bib
	pdflatex repr.tex
	bibtex repr
	pdflatex repr.tex
	pdflatex repr.tex





spectrum: spectrum.tex refs.bib
	pdflatex spectrum.tex
	bibtex spectrum
	pdflatex spectrum.tex
	pdflatex spectrum.tex
	open spectrum.pdf 


