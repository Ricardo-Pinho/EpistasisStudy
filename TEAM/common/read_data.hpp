#ifndef READ_DATA_H
#define READ_DATA_H

using namespace std;


#include "global.hpp"

int read_bi_genotypes(ifstream &geno_infile, vector<biSNP> &X);

int read_tri_genotypes(ifstream &geno_infile, vector<triSNP> &X);

int read_bi_phenotypes(ifstream &pheno_infile, biPhenotype &Y);

int read_qt_phenotypes(ifstream &pheno_infile, qtPhenotype &Y);

#endif
