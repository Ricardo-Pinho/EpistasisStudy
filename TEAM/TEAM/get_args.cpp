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

#include "global_param.hpp"

extern char* geno_fn;
extern char* pheno_fn;
extern char* pvalue_fn;
extern int MAX_NUM_OF_INDIVIDUALS;
extern int MAX_NUM_OF_SNPS;
extern int NUM_OF_PERMUTATIONS;
extern int id;
extern int n_partitions;
extern int n_threads;
extern ifstream geno_infile;
extern ifstream pheno_infile;

extern char* qvalue_fn;
extern ifstream qvalue_infile;

double qvalue_threshold;


void usage(){
	cout <<"Usage: test_all [OPTIONS]..." <<endl <<endl
         <<"File options:" <<endl
         <<"-if_geno geno_file "<<endl
         <<"-if_pheno pheno_file" <<endl
         <<"-of_pvalue pvalue_file" <<endl <<endl
         <<"Parameter options:" <<endl
         <<"-n_inds #individuals" <<endl
         <<"-n_snps #snps" <<endl
         <<"-n_perms #permutations" <<endl <<endl
//         <<"Parallel options:" <<endl
//         <<"-n_partitions #parititions" <<endl
//         <<"-id #id" <<endl
//         <<"-n_threads #threads" <<endl
         <<endl;    
    exit(1);
}


int get_args(int argc, char * argv[]) {
    if (argc<3) {
		usage();
	}	
	int arg_cur=1;
	while(arg_cur<argc){
		if (strcmp(argv[arg_cur],"-if_geno")==0){
            arg_cur++;
			if (arg_cur<argc){
				geno_fn=argv[arg_cur];
				arg_cur++;
			}
			else {
				usage();
			}
		}
        if (strcmp(argv[arg_cur],"-if_pheno")==0){
			arg_cur++;
			if (arg_cur<argc){
				pheno_fn=argv[arg_cur];
				arg_cur++;
			}
			else {
				usage();;
			}
		}
		else if (strcmp(argv[arg_cur],"-of_pvalue")==0){
			arg_cur++;
			if (arg_cur<argc){
				pvalue_fn=argv[arg_cur];
				arg_cur++;
			}
			else {
				usage();
			}
		}
        else if (strcmp(argv[arg_cur],"-n_inds")==0){
			arg_cur++;
			if (arg_cur<argc){
				MAX_NUM_OF_INDIVIDUALS=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage();
			}
		}
		else if (strcmp(argv[arg_cur],"-n_snps")==0){
			arg_cur++;
			if (arg_cur<argc){
				MAX_NUM_OF_SNPS=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage();
			}
		}		
		else if (strcmp(argv[arg_cur],"-n_perms")==0){
			arg_cur++;
			if (arg_cur<argc){
				NUM_OF_PERMUTATIONS=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage();
			}
		}
        // else if (strcmp(argv[arg_cur],"-n_partitions")==0){
		// 	arg_cur++;
		// 	if (arg_cur<argc){
		// 		n_partitions=atoi(argv[arg_cur]);
		// 		arg_cur++;
		// 	}
		// 	else {
		// 		usage();
		// 	}
		// }
        // else if (strcmp(argv[arg_cur],"-id")==0){
		// 	arg_cur++;
		// 	if (arg_cur<argc){
		// 		id=atoi(argv[arg_cur]);
		// 		arg_cur++;
		// 	}
		// 	else {
		// 		usage();
		// 	}
		// }
        // else if (strcmp(argv[arg_cur],"-n_threads")==0){
		// 	arg_cur++;
		// 	if (arg_cur<argc){
		// 		n_threads=atoi(argv[arg_cur]);
		// 		arg_cur++;
		// 	}
		// 	else {
		// 		usage();
		// 	}
		// }
        
		else {
			usage();
		}
	}

	if (geno_fn==NULL || pheno_fn==NULL) {
		usage();
		
	}	
	geno_infile.open(geno_fn);
    pheno_infile.open(pheno_fn);
    
	if(geno_infile.fail() || pheno_infile.fail()) {
		cerr <<"Cannot open file" <<endl;		
		usage();
	}
    return 0;
    
}


void usage_get_snps(){
	cout <<"Usage: get_snps [OPTIONS]..." <<endl <<endl
         <<"File options:" <<endl
         <<"-if_geno geno_file "<<endl
         <<"-if_pheno pheno_file" <<endl
         <<"-if_qvalue qvalue_file" <<endl
         <<"Parameter options:" <<endl
         <<"-n_inds #individuals" <<endl
         <<"-n_snps #snps" <<endl
         <<"-n_perms #permutations" <<endl
         <<"-qvalue qvalue_threshold" <<endl <<endl
        /*
         <<"Parallel options:" <<endl
         <<"-n_partitions #parititions" <<endl
         <<"-id #id" <<endl
         <<"-n_threads #threads" <<endl
        */
         <<endl;    
    exit(1);
}


int get_args_get_snps(int argc, char * argv[]) {
    if (argc<3) {
		usage_get_snps();
	}	
	int arg_cur=1;
	while(arg_cur<argc){
		if (strcmp(argv[arg_cur],"-if_geno")==0){
			arg_cur++;
			if (arg_cur<argc){
				geno_fn=argv[arg_cur];
				arg_cur++;
			}
			else {
				usage_get_snps();
			}
		}
        if (strcmp(argv[arg_cur],"-if_pheno")==0){
			arg_cur++;
			if (arg_cur<argc){
				pheno_fn=argv[arg_cur];
				arg_cur++;
			}
			else {
				usage_get_snps();;
			}
		}
        else if (strcmp(argv[arg_cur],"-n_inds")==0){
			arg_cur++;
			if (arg_cur<argc){
				MAX_NUM_OF_INDIVIDUALS=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage();
			}
		}
		else if (strcmp(argv[arg_cur],"-n_snps")==0){
			arg_cur++;
			if (arg_cur<argc){
				MAX_NUM_OF_SNPS=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage_get_snps();
			}
		}
        else if (strcmp(argv[arg_cur],"-n_partitions")==0){
			arg_cur++;
			if (arg_cur<argc){
				n_partitions=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage_get_snps();
			}
		}
        else if (strcmp(argv[arg_cur],"-id")==0){
			arg_cur++;
			if (arg_cur<argc){
				id=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage_get_snps();
			}
		}
        else if (strcmp(argv[arg_cur],"-n_threads")==0){
			arg_cur++;
			if (arg_cur<argc){
				n_threads=atoi(argv[arg_cur]);
				arg_cur++;
			}
			else {
				usage_get_snps();
			}
		}
        else if (strcmp(argv[arg_cur],"-if_qvalue")==0){
			arg_cur++;
			if (arg_cur<argc){
				qvalue_fn=argv[arg_cur];
				arg_cur++;
			}
			else {
				usage_get_snps();
			}
		}
        else if (strcmp(argv[arg_cur],"-qvalue")==0){
			arg_cur++;
			if (arg_cur<argc){
				qvalue_threshold=atof(argv[arg_cur]);
				arg_cur++;
                
			}
			else {
				usage_get_snps();
			}
		}
		else {
			usage_get_snps();
		}
	}

	if (geno_fn==NULL || pheno_fn==NULL || qvalue_fn==NULL) {
		usage_get_snps();
	}	
	geno_infile.open(geno_fn);
    pheno_infile.open(pheno_fn);
    qvalue_infile.open(qvalue_fn);
    
	if(geno_infile.fail() || pheno_infile.fail() || qvalue_infile.fail()) {
		cerr <<"Cannot open file" <<endl;		
		usage_get_snps();
	}
    return 0;
}
