SHELL = /bin/bash
#-------------------------------------------------------------------------------
.PHONY:	all adm bssn code
#-------------------------------------------------------------------------------
all:
	@ echo "> make adm bssn pdf ..."
	@ make adm
	@ make bssn
	@ pdflatex -halt-on-error -interaction=batchmode adm-bssn-plots &> adm-bssn-plots.texlog
	@ make veryclean
#-------------------------------------------------------------------------------
adm:
	@ echo "> make adm ..."
	@ (cd adm/code;     make; make data; make results)
#-------------------------------------------------------------------------------
bssn:
	@ echo "> make bssn ..."
	@ (cd bssn/code;    make; make data; make results)
#-------------------------------------------------------------------------------
code:
	@ echo "> make adm/code ..."
	@ (cd adm/code;     make)
	@ echo "> make bssn/code ..."
	@ (cd bssn/code;    make)
#-------------------------------------------------------------------------------
rm-dot:
	@ rm -rf .[a-z]*.lb
#-------------------------------------------------------------------------------
clean:
	@ rm -rf *.aux *.log *.out *.texlog *.synctex.gz
#-------------------------------------------------------------------------------
veryclean:
	@ make clean
#-------------------------------------------------------------------------------
pristine:
	@ make rm-dot
	@ make veryclean
	@ rm -rf adm-bssn-*.pdf
	@ (cd adm/code;     make pristine)
	@ (cd bssn/code;    make pristine)
	@ (cd support;      rm -rf obj)
#-------------------------------------------------------------------------------
github-clean:
	@ # same as "pristine" but keep the final pdf
	@ make rm-dot
	@ make veryclean
	@ (cd adm/code;     make pristine)
	@ (cd bssn/code;    make pristine)
	@ (cd support;      rm -rf obj)
#-------------------------------------------------------------------------------
github:
	@ make pristine
	@ make all
	@ make github-clean
