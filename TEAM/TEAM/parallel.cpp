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

#include "global_param.hpp"

bool get_snps_partition(__int64 n_snps, __int64 n_partitions, int partition_id, 
                       __int64 & from_pos, __int64 & to_pos, __int64 & size) // from_pos 
{
    __int64 total_pairs=n_snps*(n_snps-1)/2;
    __int64 partition_size=total_pairs/n_partitions;
    // cout <<"Total pairs: " <<total_pairs <<endl
    //      <<"Parition size: " <<partition_size <<endl;
    // cout <<"============" <<endl;
    int pre_pos=0;
    int cur_pos=0;
    __int64 cur_count=0;
    __int64 tot_count=0;
    for (int i=0 ; i<n_partitions ; ++i){
        cur_count=0;
        for (int j=n_snps-1-cur_pos ; j>=0 ; --j){
            cur_count+=j;
            cur_pos++;   
            if (cur_count>=partition_size) {                
                break;
            }         
        }
        if (cur_count<=0)
            return false; //Fail to making even partitions.
        
        //cout <<pre_pos <<"--" <<(cur_pos-1) <<"\t" <<cur_count <<endl;
        if (i==partition_id) {
            from_pos=pre_pos;
            if (cur_pos==n_snps-1)
                to_pos=cur_pos;
            else
                to_pos=cur_pos-1;
            size=cur_count;
            return true;
        }
        
        tot_count+=cur_count;
        pre_pos=cur_pos;
    }
    
    // cout <<"============="
    //      <<endl
    //      <<tot_count
    //      <<endl;
    return false;
}
