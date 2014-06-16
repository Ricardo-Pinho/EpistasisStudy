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

/*#include "global_param.hpp"*/

typedef long long __int64; //for g++

char *fn=NULL;
__int64 *bins;

const int n_bins=1000;
const double min_value=0;
const double max_value=1;

int main(int argc, char * argv[])
{
    fn=argv[1];
    
    ifstream fin(fn);
    if (fin.fail()) {
        cerr <<"ERROR: Cannot open file \"" <<fn <<"\"" <<endl;
        return 1;        
    }

    bins=new __int64[n_bins];
    memset(bins, 0, n_bins*sizeof(__int64));
    
    double col1;
    __int64 col2;

    double interval=(max_value-min_value)/n_bins;
    
    int i;

    __int64 perm_total=0;
    i=0;    
    while(fin >>col1 >>col2){
        int pos=int((col1-min_value)/interval);        
        bins[pos]+=col2;        
        ++i;
    }
    
    for (i=0 ; i<n_bins ; ++i){
        cout <<i*interval+min_value
             <<"\t" <<bins[i]
             <<endl;        
    }

    delete []bins;
}
