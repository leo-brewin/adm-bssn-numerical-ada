SHELL = /bin/sh
.PHONY:	all data merge results clean veryclean pristine
#-------------------------------------------------------------------------------
all:
	@ build.sh
data:
	@ adminitial.sh
results:
	@ admevolve.sh
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
