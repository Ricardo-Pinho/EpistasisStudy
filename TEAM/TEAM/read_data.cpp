using namespace std;

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <bitset>
#include <cstring>

#include<unistd.h>
#include<sys/types.h>
#include<stdio.h>
#include<stdlib.h>
#include<limits.h>

#include "../common/global.hpp"
#include "global_param.hpp"

extern int NUM_OF_ONES_IN_Y;
extern int MAX_NUM_OF_SNPS;

/*
int read_genotypes(ifstream &geno_infile, vector<SNP> &X)
{
    SNP temp_snp;    
    char* temp_str;
	temp_str=new char[MAX_NUM_OF_INDIVIDUALS+1];
    
	cout <<"Reading Data...";
	cout.flush();	
	for (int i = 0; i < MAX_NUM_OF_SNPS ; ++i) {
		geno_infile >>temp_str;		
		for (int j=0 ; j< MAX_NUM_OF_INDIVIDUALS ; ++j){
			switch (temp_str[j]){
			case '0': temp_snp.zero.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;
			case '1': temp_snp.one.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;				
			case '2': temp_snp.two.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;
			default: cerr <<"SNP data error!" <<endl; return 1;				
			}			
		}		
		X.push_back(temp_snp);

		// cout <<temp_snp.zero.count() <<"\t"
		// 	 <<temp_snp.one.count() <<"\t"
		// 	 <<temp_snp.two.count() <<"\t"
		// 	 <<endl;
		
		// cout <<endl <<i <<")\t" <<temp_str <<endl;
		// cout <<"-(" <<temp_snp.zero.count() <<")\t" <<temp_snp.zero <<endl
		// 	 <<"-(" <<temp_snp.one.count() <<")\t" <<temp_snp.one <<endl
		// 	 <<"-(" <<temp_snp.two.count() <<")\t" <<temp_snp.two <<endl;

		temp_snp.zero.reset();
		temp_snp.one.reset();
		temp_snp.two.reset();		
	}
	// cout <<endl <<endl;
	
	cout <<"OK" <<endl;
	cout <<X.size() <<" SNPs read." <<endl <<endl;
    return 0;    
}


Phenotype read_phenotypes(ifstream &pheno_infile) {    
    Phenotype Y_origin;
    int tmp1;
    
    for (int i=0 ; i < MAX_NUM_OF_INDIVIDUALS ; ++i){        
        pheno_infile >>tmp1;
        if (tmp1==1) {
            Y_origin.set(i);            
        }
    }
    
    //Ylist.push_back(Y_origin);
    cout <<"Original Phenotypes:" <<endl
         <<"#Cases:   \t" <<Y_origin.count() <<endl
         <<"#Controls:\t" <<MAX_NUM_OF_INDIVIDUALS-Y_origin.count() <<endl <<endl;
    NUM_OF_ONES_IN_Y=Y_origin.count();
    return Y_origin;    
}
*/
double read_qvalue_list(ifstream &qvalue_infile, double qvalue_threshold) { //return the test score threshold
    double qvalue,score;
    int count;
    cout <<"QValue Threshold:" <<qvalue_threshold <<endl;
    
    while(qvalue_infile >>qvalue >>count >>score){
        if (qvalue<=qvalue_threshold)
            return score;
    }
    return -1;
}

void generate_permutation(int num_of_permutations, int num_of_ones, list <biPhenotype>& Ylist)
{
	set<int> ones_index;	
	int r,k;	
	srand(1);
	//srand((unsigned)time( NULL ) );
	
	biPhenotype tempY(MAX_NUM_OF_INDIVIDUALS);
	
	k=0;	
	while(k < num_of_permutations) {
		for (int i = 0; i < num_of_ones;) {		
			r=rand()%MAX_NUM_OF_INDIVIDUALS;
			if (ones_index.find(r)==ones_index.end()){
				ones_index.insert(r);
				i++;				
			}		
		}
	
		for (set<int>::iterator j=ones_index.begin() ; j!=ones_index.end() ; j++){		
			//cout <<(*j) <<endl;
			tempY.set((*j));
		}
		#ifdef verbose
		cout <<tempY <<endl;
		#endif
		Ylist.push_back(tempY);
		tempY.reset();
		ones_index.clear();
		k++;		
	}
}
