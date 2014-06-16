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

char *perm_fn=NULL;
char *orig_fn=NULL;

__int64 *bins_perm;
__int64 *bins_orig;

int main(int argc, char * argv[])
{
    if (argc!=3) {
        cout <<"Usage: pvalue bins_perm_file bins_origin_file\n";        
        return 1;        
    }

    perm_fn=argv[1];
    orig_fn=argv[2];
    
    ifstream perm_in(perm_fn);

    if (perm_in.fail()) {
        cerr <<"ERROR: Cannot open file \"" <<perm_fn <<"\"" <<endl;
        return 1;        
    }

    ifstream orig_in(orig_fn);    
    if (orig_in.fail()) {
        cerr <<"ERROR: Cannot open file \"" <<orig_fn <<"\"" <<endl;
        return 1;        
    }

    bins_perm=new __int64[n_bins];
    memset(bins_perm, 0, n_bins*sizeof(__int64));

    bins_orig=new __int64[n_bins];
    memset(bins_orig, 0, n_bins*sizeof(__int64));

    
    string fn;
    double col1;
    __int64 col2;

    int i;

    __int64 perm_total=0;
    i=0;    
    while(perm_in >>col1 >>col2){
        bins_perm[i]=col2;        
        perm_total+=bins_perm[i];
        ++i;
    }
    
    for (int i=n_bins-2 ; i>=0 ; --i)        
        bins_perm[i]=bins_perm[i+1]+bins_perm[i];    

    __int64 orig_total=0;
    i=0;    
    while(orig_in >>col1 >>col2){
        bins_orig[i]=col2;            
        orig_total+=bins_orig[i];
        ++i;        
    }

    for (int i=n_bins-2 ; i>=0 ; --i)        
        bins_orig[i]=bins_orig[i+1]+bins_orig[i];    

    //cout <<bins_perm[0] <<"\t" <<bins_orig[0] <<endl;

    double p_correct;

    __int64 old_orig=bins_orig[0];
    
    for (i=0 ; i<n_bins ; ++i){
        if (bins_orig[i]==0)
            break;
        if (old_orig!=bins_orig[i]) {            
            cout <<p_correct
                 <<"\t" <<old_orig-bins_orig[i]
                 <<"\t" <<bins_perm[i-1] <<"\t" <<bins_perm[0]
                 <<"\t" <<bins_orig[i-1] <<"\t" <<bins_orig[0]
                 <<endl;
            old_orig=bins_orig[i];            
        }
        //p_correct=double(bins_orig[i])/(double(bins_perm[i]));
	p_correct=double(bins_perm[i])/double(perm_total);
    }
    cout <<p_correct
         <<"\t" <<bins_orig[i-1]
         <<"\t" <<old_orig
         <<"\t" <<bins_perm[i-1] <<"\t" <<bins_perm[0]
         <<"\t" <<bins_orig[i-1] <<"\t" <<bins_orig[0]
         <<endl;        

/*    
    __int64 orig_count=0;
    
    for (i=0 ; i<n_bins ; ++i){
        if (bins_orig[i]==0)
            p_correct=0;
        else 
            p_correct=double(bins_perm[i])/(perm_total)/(double(bins_orig[i])/(orig_total));

        if (orig_count>=orig_total)
            break;
        orig_count+=bins_orig[i];
        cout <<p_correct <<"\t" <<bins_orig[i]
             <<endl;
    }
*/    
    delete []bins_orig;
    delete []bins_perm;    
}
