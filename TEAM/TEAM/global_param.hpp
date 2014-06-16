#ifndef GLOBAL_PAR_H
#define GLOBAL_PAR_H


/* Human data */
//const int MAX_NUM_OF_INDIVIDUALS=500;
//const int MAX_NUM_OF_INDIVIDUALS=2520;
//const int MAX_NUM_OF_INDIVIDUALS=3073;


/*
struct SNP {
    bitset<MAX_NUM_OF_INDIVIDUALS> zero,one,two;    
};

typedef bitset<MAX_NUM_OF_INDIVIDUALS> Phenotype;
typedef bitset<MAX_NUM_OF_INDIVIDUALS> Intermedia;
*/
typedef map < pair<int,int> , vector<int> >	Data_map;

typedef long long __int64; //for g++

struct SNPs_Diff
{
	set <short> zero2one,zero2two,
		one2zero,one2two,
		two2zero,two2one;
};

#endif
