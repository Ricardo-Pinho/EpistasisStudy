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


#include<unistd.h>
#include<sys/types.h>
#include<stdio.h>
#include<stdlib.h>
#include<limits.h>

#include "global_param.hpp"
#include "mst.hpp"

extern int MAX_NUM_OF_SNPS;

bool operator<(const edge& a, const edge& b){
	if ((a.weight) == (b.weight))
		if ((a.p1)==(b.p1)) return (a.p2)>(b.p2);	
		else return (a.p1)>(b.p1);	
	else return (a.weight) > (b.weight);
}

bool operator>(const edge& a, const edge& b){
	if ((a.weight) == (b.weight))
		if ((a.p1)==(b.p1)) return (a.p2)>(b.p2);	
		else return (a.p1)>(b.p1);	
	else return (a.weight) < (b.weight);
}

//Calculate distance for a SNP-pair.
int calculate_distance(triSNP &s1, triSNP &s2){
	triSNP tempSNP;
	tempSNP.zero=s1.zero & s2.zero;
	tempSNP.one=s1.one & s2.one;
	tempSNP.two=s1.two & s2.two;
	
	return MAX_NUM_OF_INDIVIDUALS-tempSNP.zero.count()-tempSNP.one.count()-tempSNP.two.count();
}

//Calculate distance matrix for a subset of SNPs.
void calculate_dist_matrix(short dist_matrix[MAX_SUBMST_SIZE][MAX_SUBMST_SIZE], vector<triSNP>::iterator X_begin, vector<triSNP>::iterator X_end) {
	vector<triSNP>::iterator i,j;
	int ii,jj,w;	
	triSNP tempSNP;	
	for (i = X_begin, ii=0; i != X_end; i++,ii++){
		for (j = i, j++, jj=ii+1; j != X_end ; j++,jj++){
			tempSNP.zero=(*i).zero & (*j).zero;
			tempSNP.one=(*i).one & (*j).one;
			tempSNP.two=(*i).two & (*j).two;
			
			dist_matrix[ii][jj]=MAX_NUM_OF_INDIVIDUALS-(tempSNP.zero.count()+tempSNP.one.count()+tempSNP.two.count());
			dist_matrix[jj][ii]=dist_matrix[ii][jj];			
		}		
	}
}


//Build Minimum Spanning Tree for a subset of SNPs
int build_sub_mst(vector<triSNP>::iterator X_begin, vector<triSNP>::iterator X_end, list<edge> & mst, int submst_size, int offset) {	
	int new_vertex=0;
	int mst_total_weight=0;
	vector<edge> edges;
	set<short> mst_vertex;
	short dist_matrix[MAX_SUBMST_SIZE][MAX_SUBMST_SIZE]={0};	
		
	/* Calculating Distance*/
	
	clock_t t_tmp = clock();
	calculate_dist_matrix(dist_matrix,X_begin,X_end);
	clock_t t_calc_dist = clock();

	
	// for (vector<edge>::iterator j=edges.begin() ; j!=edges.end() ; j++){
	// 	cout <<"(" <<j->p1 <<"," <<j->p2 <<"):" <<j->weight <<endl;		
	// }
	
	// cout <<MAX_NUM_OF_SNPS <<" vertexes." <<endl;	
	// cout <<edges.size() <<" edges." <<endl;	
	// cout <<"Done!" <<endl <<"Time Elapsed for Calculating Distance:" <<(t_calc_dist-t_tmp)*1.0/CLOCKS_PER_SEC <<" s" <<endl <<endl;

	
/*Prim's Algorithm*/	
	mst_vertex.insert(new_vertex);	
	make_heap(edges.begin(),edges.end());
	
	while (mst_vertex.size()<submst_size){
		//Push edges and maintain the heap
		for (int k = 0 ; k<submst_size ; ++k){
			if (mst_vertex.find(k)==mst_vertex.end() && new_vertex!=k){
				//Vertex not included in current MST

				//Make sure p1<p2.
				if (new_vertex<k)
					edges.push_back(edge(new_vertex,k,dist_matrix[new_vertex][k]));
				else
					edges.push_back(edge(k,new_vertex,dist_matrix[new_vertex][k]));
				push_heap(edges.begin(),edges.end());
			}
		}
		
		while(edges.size()>0){			
			pop_heap(edges.begin(),edges.end());
			edge tempEdge(edges.back());		
			edges.pop_back();

			if (mst_vertex.find(tempEdge.p1)==mst_vertex.end()){
				new_vertex=tempEdge.p1;				
				tempEdge.p1=tempEdge.p1+offset;
				tempEdge.p2=tempEdge.p2+offset;				
				mst.push_back(tempEdge);
				mst_vertex.insert(new_vertex);
				mst_total_weight+=tempEdge.weight;				
				break;				
			}
			else if (mst_vertex.find(tempEdge.p2)==mst_vertex.end()) {
				new_vertex=tempEdge.p2;
				tempEdge.p1=tempEdge.p1+offset;
				tempEdge.p2=tempEdge.p2+offset;				
				mst.push_back(tempEdge);
				mst_vertex.insert(new_vertex);
				mst_total_weight+=tempEdge.weight;
				break;
			}
		}
	}
/*
		 for (list<edge>::iterator l=mst.begin() ; l!=mst.end() ; l++){
		 	cout <<"(" <<l->p1 <<"," <<l->p2 <<"):" <<l->weight <<"; ";
		 }
		 cout <<endl;
*/
	return mst_total_weight;	
}


//Divide SNPs into subsets, build MST for each subsets, and connect them to get a large tree,
int build_mst(vector<triSNP>& X, list<edge> & mst) {
	int offset;
	int submst_size;
	int mst_total_weight=0;
	vector<triSNP>::iterator X_begin,X_end;
    
    srand(2);
	X_begin=X.begin();	
	for (offset=0 ; offset<MAX_NUM_OF_SNPS ; offset+=MAX_SUBMST_SIZE, X_begin+=MAX_SUBMST_SIZE){
		if ( (offset+MAX_SUBMST_SIZE)<=MAX_NUM_OF_SNPS ){
			submst_size=MAX_SUBMST_SIZE;
			X_end=X_begin;
			X_end+=MAX_SUBMST_SIZE;			
		}
		else {
			submst_size=MAX_NUM_OF_SNPS-offset;
			X_end=X.end();			
		}
		//Build a sub mst
		mst_total_weight+=build_sub_mst(X_begin,X_end,mst,submst_size,offset);
		
		if (offset==0) continue;
		else {
			//Connect the sub mst to the overall mst.
			int pseudo_min_dist=0x7fffffff;
			
			int p1,p2;			
			for (int j = offset ; j < offset+submst_size ; ++j){
				for (int i = 0 ; i<CONNECTION_TRIALS ; ++i){
					int r=rand()%offset;
					int dist=calculate_distance(X[r],X[j]);
					if (pseudo_min_dist>dist){
						p1=r;
						p2=j;
						pseudo_min_dist=dist;						
					}					
				}
			}
			mst.push_back(edge(p1,p2,pseudo_min_dist));
			mst_total_weight+=pseudo_min_dist;
		}		
	}
	return mst_total_weight;	
}

//Print Function
void print_mst(list<edge> & mst){
	for (list<edge>::iterator l=mst.begin() ; l!=mst.end() ; l++){
		cout <<"(" <<l->p1 <<"," <<l->p2 <<"):" <<l->weight <<"; ";
	}
	cout <<endl;
}

//Generate a linear tree regardless its weight.
int build_linear_tree(vector<triSNP>& X, list<edge> & mst) {
	int mst_total_weight=0;	
	for (int i=0 ; i<(X.size()-1) ; ++i){
		int weight=calculate_distance(X[i],X[i+1]);
		mst.push_back(edge(i,i+1,weight));
		mst_total_weight+=weight;		
	}
	
	return mst_total_weight;	
}
