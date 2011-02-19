#!/usr/bin/env python
"""
Functions to use Funmap package to analyze mouse activity data.
"""

import csv
import numpy as np
from rpy2 import robjects as robj

def derive_cM(input_fl="../data/N1_GenoTypes.csv", output_fl="N1_Geno_cM"):
    """
    Compute the cM from basepair.
    """
    R = robj.r
    R("""
    library(qtl)
    cr <- read.cross(file='%s', format='csv', alleles=c('A', 'H'))
    #cr <- rescalemap(cr)
    mp <- est.map(cr, offset=0)
    cr <- replace.map(cr,mp)
    cr <- subset(cr,ind=1:89)
    cr <- subset(cr,chr=1:19)
    write.cross(cr, format='csv', filestem='%s')
    """ % (input_fl, output_fl))

def make_geno_files(input_fl="N1_Geno_cM.csv", output_marker="N1_marker.csv", output_geno="N1_geno.csv"):
    """
    Extract genotypes, marker info. from the input file and write to the output files.
    """
    with open(input_fl) as fl:
        stream = csv.reader(fl, delimiter=',')
        ## extract marker info
        names = next( stream)[2:]
        chrom =  next(stream)[2:]
        basepairs =  [float(x) for x in next(stream)[2:] if len(x)>0]
        assert len(names) == len(chrom) == len(basepairs), "Marker info are not of same length."
        ## extract genotype info
        geno_lst = [[int(x[0])]+list(encode(x[2:])) for x in stream if len(x)>0]
                        
    ## write output files
    with open(output_marker, 'w') as marker:
        dest = csv.writer(marker, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        lst = zip(names, names,basepairs, chrom, ['G'+str(x) for x in chrom])
        dest.writerow( ['Id', 'Marker', 'Distance', 'chromosome(group) index', 'chromosome(group)name'] ) # write the header
        dest.writerows(lst)
    with open(output_geno, 'w') as geno:
        dest = csv.writer(geno, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        dest.writerow( ['ID'] + [str(x) for x in range(len(names))] ) # write the header
        dest.writerows( geno_lst)

def encode(seq):
    for x in seq:
        if x == 'AA':
            yield 0
        elif x == 'AH':
            yield 1
        elif x == '-':
            yield -1
        else:
            raise Exception("Unknown genotype")

def get_times(flnm="../data/N1_BinNumToCT_BW6min.csv"):
    with open(flnm) as fl:
        dt = csv.reader(fl, delimiter=',', quotechar='"')
        assert len( next(dt) ) == 222, "Incorrect bin numbers" # useless header
        start = np.array(next(dt), dtype=float)
        end = np.array(next(dt), dtype=float)
        assert start.shape == end.shape == (222,), "Start or end time-points are not correct shape."
    tmp = (start+end)/2.
    return tmp

def make_pheno_file(input_fl="../data/N1_CTvsProb_BW6min.csv", output_fl="N1_pheno.csv"):
    """
    Convert to phenotype file format for Funmap.
    """
    tp = get_times()
    with open(input_fl) as fl:
        istream = csv.reader(fl, delimiter=',')
        next(istream)                    # useless header
        names = list(); pheno = list
        with open(output_fl, 'w') as ofl:
            ostream = csv.writer(ofl, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
            ostream.writerow(['ID']+get_times().tolist())
            [ostream.writerow([int(x[0])]+[float(y) for y in x[1:]]) for x in istream if len(x)>0] # write out non-empty lines

if __name__ == '__main__':
    derive_cM()
    make_geno_files()
    make_pheno_file()
