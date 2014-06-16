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


typedef long long __int64; //for g++

char *bin_fn=NULL;

const int n_bins=10000000;

double *bins_score;
__int64 *bins_perm;
__int64 *bins_orig;


int main(int argc, char * argv[])
{
    if (argc!=2) {
        cout <<"Usage: pvalue collected_bins\n";        
        return 1;        
    }

    bin_fn=argv[1];
    
    ifstream bin_in(bin_fn);

    if (bin_in.fail()) {
        cerr <<"ERROR: Cannot open file \"" <<bin_fn <<"\"" <<endl;
        return 1;        
    }

    bins_score=new double[n_bins];
    memset(bins_score, 0, n_bins*sizeof(double));
    
    bins_perm=new __int64[n_bins];
    memset(bins_perm, 0, n_bins*sizeof(__int64));

    bins_orig=new __int64[n_bins];
    memset(bins_orig, 0, n_bins*sizeof(__int64));
    
    string fn;
    double col1;
    __int64 col2;
    __int64 col3;

    int i;
    
    i=0;    
    while(bin_in >>col1 >>col2 >>col3){
        bins_score[i]=col1;
        bins_perm[i]=col2;            
        bins_orig[i]=col3;
        ++i;
    }

    
    for (int i=n_bins-2 ; i>=0 ; --i){        
        bins_perm[i]=bins_perm[i+1]+bins_perm[i]; 
        bins_orig[i]=bins_orig[i+1]+bins_orig[i];    
    }    

    //cout <<bins_perm[0] <<"\t" <<bins_orig[0] <<endl;

    //cout <<"OK" <<endl;

    double p_correct;

    __int64 old_orig=bins_orig[0];
    
    for (i=0 ; i<n_bins-1 ; ++i){
        if (bins_orig[i]==0)
            break;
        if ( (bins_orig[i]-bins_orig[i+1])>0 ) {            
            cout <<p_correct
                //<<"\t" <<bins_perm[i-1] //<<"\t" <<bins_perm[0]
                 <<"\t" <<bins_orig[i-1]-bins_orig[i+1] //<<"\t" <<bins_orig[0]
                 <<"\t" <<bins_score[i]
                 <<endl;
        }
        p_correct=double(bins_perm[i])/(bins_perm[0])/(double(bins_orig[i])/(bins_orig[0]));
    }
    if (bins_orig[n_bins-1]>0) {
        cout <<p_correct 
            //<<"\t" <<bins_perm[i-1] //<<"\t" <<bins_perm[0]
             <<"\t" <<bins_orig[n_bins-1] //<<"\t" <<bins_orig[0]
             <<"\t" <<bins_score[n_bins-1]
             <<endl;
    }

    //cout <<bin_interval <<endl;
    delete []bins_orig;
    delete []bins_perm;    
}
