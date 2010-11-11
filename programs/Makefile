.PHONY: clean purge programfiles all
# group files by type
Rfiles = R/fr.R R/logistic.R 
pyfiles :=  $(wildcard py/*.py) $(wildcard py/*.pyx) py/Makefile py/README.txt
datafiles = data/N1_GenoTypes.csv data/N1_CTvsProb_BW6min.csv

all: xiongSupplement.tgz

programfiles: functionalMappingPrograms.tgz

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

functionalMappingPrograms.tgz: $(Rfiles) $(pyfiles)
        # create directories
	mkdir -p functionalMappingPrograms/R
	mkdir functionalMappingPrograms/py
# copy files into the respective directories
	cp $(Rfiles) functionalMappingPrograms/R
	cp $(pyfiles) functionalMappingPrograms/py
# tarzip it
	tar -zcf functionalMappingPrograms.tgz functionalMappingPrograms


clean:
	rm -rf xiongSupplement
	rm -rf functionalMappingPrograms

purge:
	rm -rf xiongSupplement
	rm -f xiongSupplement.tgz
	rm -rf functionalMappingPrograms
	rm -f functionalMappingPrograms.tgz
