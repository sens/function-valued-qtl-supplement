# gzipped tar file dependencies
xiongSupplement.tgz: R/fr.R R/logistic.R R/steps1.R py/*.{py,pyx} py/Makefile\
	data/N1_GenoTypes.csv data/N1_CTvsProb_BW6min.csv
	# create directories
	mkdir -p xiongSupplement/R
	mkdir xiongSupplement/py
	mkdir xiongSupplement/data
# copy files into the respective directories
	cp R/fr.R R/logistic.R xiongSupplement/R
	cp data/N1_GenoTypes.csv xiongSupplement/data
	cp data/N1_CTvsProb_BW6min.csv xiongSupplement/data
	cp py/*.{py,pyx} py/Makefile xiongSupplement/py
# tarzip it
	tar -zcf xiongSupplement.tgz xiongSupplement

clean:
	rm -rf xiongSupplement
