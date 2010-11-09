.PHONY: clean purge
# group files by type
Rfiles = R/fr.R R/logistic.R 
pyfiles :=  $(wildcard py/*.py) $(wildcard py/*.pyx) py/Makefile
datafiles = data/N1_GenoTypes.csv data/N1_CTvsProb_BW6min.csv

# gzipped tar file dependencies
xiongSupplement.tgz: $(Rfiles) $(pyfiles) $(datafiles)
	# create directories
	mkdir -p xiongSupplement/R
	mkdir xiongSupplement/py
	mkdir xiongSupplement/data
# copy files into the respective directories
	cp $(Rfiles) xiongSupplement/R
	cp $(datafiles) xiongSupplement/data
	cp $(pyfiles) xiongSupplement/py
# tarzip it
	tar -zcf xiongSupplement.tgz xiongSupplement

clean:
	rm -rf xiongSupplement
purge:
	rm -rf xiongSupplement
	rm xiongSupplement.tgz
