xiongSupplement.tgz: R/fr.R R/logistic.R R/steps1.R
	mkdir xiongSupplement
	mkdir xiongSupplement/R
	mkdir xiongSupplement/py
	mkdir xiongSupplement/data
	cp R/fr.R R/logistic.R xiongSupplement/R
	cp data/N1_GenoTypes.csv xiongSupplement/data
	cp data/N1_CTvsProb_BW6min.csv xiongSupplement/data
	tar zcf xiongSupplement.tgz xiongSupplement

clean:
	rm -r xiongSupplement
