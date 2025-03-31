SHELL = /bin/sh
.PHONY:	all data results bin obj clean veryclean pristine
#-------------------------------------------------------------------------------
all:
	@ build.sh
data:
	@ bssninitial.sh
results:
	@ bssnevolve.sh
#-------------------------------------------------------------------------------
clean:
	@ rm -rf obj
veryclean:
	@ make clean
pristine:
	@ make veryclean
	@ rm -rf *.log
	@ rm -rf bin
	@ rm -rf data
	@ rm -rf results/xy results/xz results/yz
