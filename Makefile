

out: spectrum
	open spectrum.pdf 

spectrum: spectrum.tex refs.bib
	pdflatex spectrum.tex
	bibtex spectrum
	pdflatex spectrum.tex
	pdflatex spectrum.tex


