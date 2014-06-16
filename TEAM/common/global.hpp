#ifndef GLOBAL_H
#define GLOBAL_H

#include <boost/dynamic_bitset.hpp>

extern int MAX_NUM_OF_INDIVIDUALS;
extern int MAX_NUM_OF_SNPS;
extern int NUM_OF_PERMUTATIONS;

typedef boost::dynamic_bitset<> biSNP;


typedef struct {
    boost::dynamic_bitset<> zero,one,two;
} triSNP;


typedef boost::dynamic_bitset<> biPhenotype;

typedef boost::dynamic_bitset<> Intermedia;

typedef std::vector<double> qtPhenotype;

typedef long long __int64;

#endif
