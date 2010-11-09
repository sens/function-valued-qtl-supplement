# gzipped tar file dependencies
xiongSupplement.tgz: R/fr.R R/logistic.R \
	data/N1_GenoTypes.csv data/N1_CTvsProb_BW6min.csv
# make directory
	mkdir xiongSupplement
# make subdirectories
	mkdir xiongSupplement/R
	mkdir xiongSupplement/py
	mkdir xiongSupplement/data
# copy files into the respective directories
	cp R/fr.R R/logistic.R xiongSupplement/R
	cp data/N1_GenoTypes.csv xiongSupplement/data
	cp data/N1_CTvsProb_BW6min.csv xiongSupplement/data
# tarzip it
	tar zcf xiongSupplement.tgz xiongSupplement

clean:
	rm -r xiongSupplement
