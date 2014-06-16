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

char *collection=NULL;
double *bin_values;
__int64 *bins;

int main(int argc, char * argv[])
{
    if (argc!=2) {
        cout <<"Usage: collector collection_file\n";        
        return 1;        
    }

    collection=argv[1];    
    ifstream collect_infile(collection);

    if (collect_infile.fail()) {
        cerr <<"ERROR: Cannot open file \"" <<collection <<"\"" <<endl;
        return 1;        
    }

    bin_values=new double[n_bins];
    memset(bin_values, 0, n_bins*sizeof(double));
    bins=new __int64[n_bins];
    memset(bins, 0, n_bins*sizeof(__int64));

    string fn;
    double col1;
    __int64 col2;    
    while (collect_infile>>fn){
        int i=0;        
        ifstream infile(fn.c_str());
        while(infile >>col1 >>col2){
            //cout <<col1 <<"\t" <<col2 <<endl;
            bin_values[i]=col1;
            bins[i]+=col2;            
            ++i;            
        }        
    }
    for (int i=0 ; i<n_bins ; ++i){
        cout <<bin_values[i] <<"\t" <<bins[i] <<endl;
    }    
    
    delete []bin_values;
    delete []bins;    
}
