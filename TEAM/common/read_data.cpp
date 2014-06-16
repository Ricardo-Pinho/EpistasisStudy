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
#include <fstream>


#include "global.hpp"

#define verbose
extern int MAX_NUM_OF_SNPS;
extern int MAX_NUM_OF_INDIVIDUALS;


int read_bi_genotypes(ifstream &geno_infile, vector<biSNP> &X) {
    biSNP temp_snp(MAX_NUM_OF_INDIVIDUALS);
    string temp_str;

#ifdef verbose    
	cout <<"Reading Data...";
	cout.flush();
#endif
    
	for (int i = 0; i < MAX_NUM_OF_SNPS ; ++i) {
		geno_infile >>temp_str;        
        
        if (temp_str.size()!=MAX_NUM_OF_INDIVIDUALS) {            
            cerr <<"Error (Line " <<i+1 <<"): The number of Individuals does not match the preset parameter!" <<endl;
            exit(1);            
            return -1;
        }
        
        // cout <<temp_str <<endl;
		for (int j=0 ; j< MAX_NUM_OF_INDIVIDUALS ; ++j){            
			switch (temp_str[j]){
			    case '0': break;
			    case '1': temp_snp.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;
//                case '1': temp_snp.set(j); break;
                default: cerr <<"Error (Line " <<i+1 <<"): Unknown charater found in the SNP data!" <<endl; return -1;
			}            
		}
        
        if (temp_snp.count()>(MAX_NUM_OF_INDIVIDUALS-temp_snp.count()))
			temp_snp.flip(); //Make sure 1 is the minor allele

		X.push_back(temp_snp);

        // cout <<temp_snp <<endl;        
		
		// cout <<endl <<i <<")\t" <<temp_str <<endl;
		// cout <<"-(" <<temp_snp.zero.count() <<")\t" <<temp_snp.zero <<endl
		// 	 <<"-(" <<temp_snp.one.count() <<")\t" <<temp_snp.one <<endl
		// 	 <<"-(" <<temp_snp.two.count() <<")\t" <<temp_snp.two <<endl;

		temp_snp.reset();		
	}    
	// cout <<endl <<endl;
    
#ifdef verbose	
	cout <<"OK" <<endl;
	cout <<X.size() <<" SNPs read." <<endl <<endl;
#endif
    
    if (MAX_NUM_OF_SNPS!=X.size()){
        cerr <<"Error: The number of SNPs does not match the preset parameter!" <<endl;
        return -1;
    }
    return X.size();
}

int read_tri_genotypes(ifstream &geno_infile, vector<triSNP> &X) {    
    triSNP temp_snp;
    temp_snp.zero.resize((MAX_NUM_OF_INDIVIDUALS));
    temp_snp.one.resize((MAX_NUM_OF_INDIVIDUALS));
    temp_snp.two.resize(MAX_NUM_OF_INDIVIDUALS);
    
    string temp_str;

#ifdef verbose    
	cout <<"Reading Data...";
	cout.flush();
#endif
    
	for (int i = 0; i < MAX_NUM_OF_SNPS ; ++i) {
		geno_infile >>temp_str;
        if (temp_str.size()!=MAX_NUM_OF_INDIVIDUALS) {            
            cerr <<"Error (Line " <<i+1 <<"): The number of Individuals does not match the preset parameter!" <<endl;
            exit(1);            
            return -1;
        }
		for (int j=0 ; j< MAX_NUM_OF_INDIVIDUALS ; ++j){
			switch (temp_str[j]){
			case '0': temp_snp.zero.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;
			case '1': temp_snp.one.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;				
			case '2': temp_snp.two.set(MAX_NUM_OF_INDIVIDUALS-1-j); break;
			default:  cerr <<"Error: Unknown charater found in the SNP data!" <<endl; return -1;
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
#ifdef verbose	
	cout <<"OK" <<endl;
	cout <<X.size() <<" SNPs read." <<endl <<endl;
#endif

    if (MAX_NUM_OF_SNPS!=X.size()){
        cerr <<"Error: The number of SNPs does not match the preset parameter!" <<endl;
        return -1;
    }
    return X.size();
}

int read_bi_phenotypes(ifstream &pheno_infile, biPhenotype &Y) {
    int num_of_ones;
    
    char* temp_str;
	temp_str=new char[MAX_NUM_OF_INDIVIDUALS+1];
    pheno_infile >>temp_str;
    //cout <<temp_str <<endl;
    
    
    for (int i=0 ; i < MAX_NUM_OF_INDIVIDUALS ; ++i)
        if (temp_str[i]=='1') Y.set(MAX_NUM_OF_INDIVIDUALS-1-i);
//        if (temp_str[i]=='1') Y.set(i);
    

#ifdef verbose    
    cout <<"Original Phenotypes:" <<endl
         <<"#Cases:   \t" <<Y.count() <<endl
         <<"#Controls:\t" <<MAX_NUM_OF_INDIVIDUALS-Y.count() <<endl <<endl;
#endif

    num_of_ones=Y.count(); //Usually, the ones in Y stand for cases.
    return num_of_ones;
}

int read_qt_phenotypes(ifstream &pheno_infile, qtPhenotype &Y) {
    double temp_float;
    
    for (int i=0 ; i < MAX_NUM_OF_INDIVIDUALS ; ++i){
        pheno_infile >>temp_float;
        Y.push_back(temp_float);
    }
    
#ifdef verbose    
    cout <<"Original Phenotypes:" <<endl;
    for (int i=0 ; i<Y.size() ; ++i)
        cout <<Y[i] <<"\t";
    cout <<endl <<endl;
#endif
    
    return 0;    
}
