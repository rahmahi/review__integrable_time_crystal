NAME=main
VIEWPDF=evince
BIBTEXSRCS=main.bib
USE_PDFLATEX=1
CLEAN_FILES+=main.synctex.gz mainNotes.bib main.fdb_latexmk main.tar.gz *.bak *.tmp main.fls main.tdo
EXTRA_DIST+=*.pdf iop* 

include /usr/share/latex-mk/latex.gmk