#ifndef MST_H
#define MST_H

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

#include "../common/global.hpp"
#include "global_param.hpp"

class edge{
public:
    int p1,p2;
    int weight;
    
    edge(int p1,int p2,int weight)
    {
	this->p1=p1;
	this->p2=p2;
	this->weight=weight;			
    }	
};

const int MAX_SUBMST_SIZE=1000;

const int CONNECTION_TRIALS=10;

bool operator<(const edge& a, const edge& b);

bool operator>(const edge& a, const edge& b);

int calculate_distance(triSNP &s1, triSNP &s2);

void calculate_dist_matrix(short dist_matrix[MAX_SUBMST_SIZE][MAX_SUBMST_SIZE], vector<triSNP>::iterator X_begin, vector<triSNP>::iterator X_end);

int build_sub_mst(vector<triSNP>::iterator X_begin, vector<triSNP>::iterator X_end, list<edge> & mst, int submst_size, int offset);

int build_mst(vector<triSNP>& X, list<edge> & mst);
int build_linear_tree(vector<triSNP>& X, list<edge> & mst);

void print_mst(list<edge> & mst);

#endif
