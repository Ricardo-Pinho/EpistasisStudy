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
#include "../common/read_data.hpp"
#include "global_param.hpp"
#include "mst.hpp"
#include "global_headers.hpp"

//#define verbose

int n_partitions=1;
int id=0;
int n_threads=1;

const int n_bins=10000000;
double bin_interval;

int MAX_NUM_OF_INDIVIDUALS=32;
int MAX_NUM_OF_SNPS=100;
int NUM_OF_PERMUTATIONS=10;

int NUM_OF_ONES_IN_Y;

char *geno_fn=NULL;
char *pheno_fn=NULL;
char *pvalue_fn=NULL;
char *qvalue_fn=NULL;

ifstream geno_infile;
ifstream pheno_infile;
ifstream qvalue_infile;

//const int n_bins=1000;  //for test;

map< pair<int,int>, SNPs_Diff> snps_diff;

list<edge> mst;

vector<triSNP> X;
list <biPhenotype> Ylist;


double **look_up_n0;
double **look_up_n1;

short **XiY_B;
short **XiY_E;

void calculate_snps_diff(vector<triSNP>& X, list<edge> & mst, map< pair<int,int>, SNPs_Diff> & snps_diff){
    //SNPs_Diff dst_diff;
    int k=0;
    __int64 mst_size=mst.size();
    __int64 progress,o_progress=0;    
    
	for (list<edge>::iterator iter=mst.begin() ; iter!=mst.end() ; iter++){
		SNPs_Diff temp_diff;
		Intermedia temp;		
		int x1=iter->p1;
		int x2=iter->p2;

		//zero2one
		temp=X[x1].zero & X[x2].one;
		for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS ; ++i){			
			if (temp.test(i)){
				temp_diff.zero2one.insert(i);				
			}
		}
		//zero2two
		temp=X[x1].zero & X[x2].two;
		for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS ; ++i){			
			if (temp.test(i)){
				temp_diff.zero2two.insert(i);				
			}
		}
		
		//one2zero
		temp=X[x1].one & X[x2].zero;
		for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS ; ++i){			
			if (temp.test(i)){
				temp_diff.one2zero.insert(i);				
			}
		}
		//one2two
		temp=X[x1].one & X[x2].two;
		for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS ; ++i){			
			if (temp.test(i)){
				temp_diff.one2two.insert(i);				
			}
		}

        //two2zero
		temp=X[x1].two & X[x2].zero;
		for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS ; ++i){			
			if (temp.test(i)){
				temp_diff.two2zero.insert(i);				
			}
		}
		//two2one
		temp=X[x1].two & X[x2].one;
		for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS ; ++i){			
			if (temp.test(i)){
				temp_diff.two2one.insert(i);				
			}
		}		
		pair<int,int> point_pair=make_pair(x1,x2);
		snps_diff.insert(make_pair(point_pair,temp_diff));
        // if (x1==4 &&x2==5)
        //     dst_diff=temp_diff;
        k++;

        progress=k*1000/mst.size();

        if ((progress-o_progress)>=1){
            o_progress=progress;            
            cout <<setw(5) <<progress/10.0 <<"%" <<"\r";
            cout.flush();            
        }        
    }
    cout <<endl;
}

int calculate_test_scores(Data_map & data_map,
                          vector<short> &a2_y,vector<short> &a3_y,
                          vector<short> &b2_y,vector<short> &b3_y,
                          vector<short> &e2_y,vector<short> &e3_y,
                          int s, int t, int r, int p, int q, int o, int g, int h, int l,
                          int **abcdef,__int64 *bins, __int64 *org_bins, double & chi_alpha, double *chi_alpha_array){

    int ind=0;    
    double real_chi1=chi_square(a2_y[ind],a3_y[ind],abcdef[0][ind],abcdef[2][ind],s,p,g);
    double real_chi2=chi_square(b2_y[ind],b3_y[ind],abcdef[1][ind],abcdef[3][ind],t,q,h);
    double real_chi3=chi_square(e2_y[ind],e3_y[ind],abcdef[4][ind],abcdef[5][ind],r,o,l);
    double result=real_chi1+real_chi2+real_chi3;
    ++org_bins[int(result/bin_interval)];
    if (chi_alpha_array[ind]<result) {        
        chi_alpha_array[ind]=result;
        if (chi_alpha<chi_alpha_array[ind]){
            chi_alpha=chi_alpha_array[ind];
        }
    }
    
	for (int ind=1 ; ind<NUM_OF_PERMUTATIONS ; ++ind){
        //cout <<ind <<":" <<endl;        
		//double real_chi1=chi_square(a2_y[l],abcd[l][0],abcd[l][2],p,NUM_OF_ONES_IN_Y);
		//double real_chi2=chi_square(b2_y[l],abcd[l][1],abcd[l][3],q,NUM_OF_ONES_IN_Y);
        
        double real_chi1=chi_square(a2_y[ind],a3_y[ind],abcdef[0][ind],abcdef[2][ind],s,p,g);
        double real_chi2=chi_square(b2_y[ind],b3_y[ind],abcdef[1][ind],abcdef[3][ind],t,q,h);
        double real_chi3=chi_square(e2_y[ind],e3_y[ind],abcdef[4][ind],abcdef[5][ind],r,o,l);
                
		double result=real_chi1+real_chi2+real_chi3;
        ++bins[int(result/bin_interval)];        
		//output <<result <<"\n";
        if (chi_alpha_array[ind]<result){
            chi_alpha_array[ind]=result;
            if (chi_alpha<chi_alpha_array[ind]){
                chi_alpha=chi_alpha_array[ind];
            }
        }
        
	}	
    //cout <<endl;    
    return 0;
}

void pvalue(ofstream &pvalue_output, __int64 *bins_perm, __int64 *bins_orig){
    for (int i=n_bins-2 ; i>=0 ; --i){        
        bins_perm[i]=bins_perm[i+1]+bins_perm[i]; 
        bins_orig[i]=bins_orig[i+1]+bins_orig[i];    
    }
    double p_correct;

    __int64 old_orig=bins_orig[0];
    
    for (int i=0 ; i<n_bins-1 ; ++i){
        if (bins_orig[i]==0)
            break;
        if ( (bins_orig[i]-bins_orig[i+1])>0 ) {            
            pvalue_output <<p_correct
                        //<<"\t" <<bins_perm[i-1] //<<"\t" <<bins_perm[0]
                          <<"\t" <<bins_orig[i-1]-bins_orig[i+1] //<<"\t" <<bins_orig[0]
                          <<"\t" <<i*bin_interval            
                          <<endl;
        }
        p_correct=double(bins_perm[i])/(bins_perm[0])/(double(bins_orig[i])/(bins_orig[0]));
    }
    if (bins_orig[n_bins-1]>0) {
        pvalue_output <<p_correct 
                    //<<"\t" <<bins_perm[i-1] //<<"\t" <<bins_perm[0]
                      <<"\t" <<bins_orig[n_bins-1] //<<"\t" <<bins_orig[0]
                      <<"\t" <<(n_bins-1)*bin_interval
                      <<endl;
    }
}
    
//Use iteration instead of recursion to traverse.
void traverse_mst(int root, int new_root, vector< set<edge> > &vertex_edges,
				  map< pair<int,int>, SNPs_Diff> & snps_diff,
				  Data_map &data_map,
				  vector <triSNP> &X,
				  short **XiY_B, short **XiY_E,
				  vector<short> &a2_y, vector<short> &a3_y, vector<short> &b2_y,  vector<short> &b3_y,  vector<short> &e2_y,  vector<short> &e3_y,
				  vector <vector<short> > &A_tran, vector<vector<short> > &B_tran, vector<vector<short> > &E_tran,
				  map< pair<int,int>, int* > &plus_map,
				  int p, int q, int o, int g, int h, int l,int **abcdef,double &chi_alpha, double *chi_alpha_array,
                  __int64 &real_weight,__int64 &ideal_weight,clock_t &tt_traverse_tree,clock_t &tt_update_old,clock_t &tt_cal_bounds,__int64 &old_update_times, __int64 *bins, __int64 *org_bins){
	
	vector <int> check_points;
	vector <int> cp_p,cp_q,cp_o,cp_g,cp_h,cp_l;
	
	vector < vector<short> > cp_a2_y, cp_a3_y, cp_b2_y, cp_b3_y, cp_e2_y, cp_e3_y;

	vector <short> current_a2_y(NUM_OF_PERMUTATIONS), current_a3_y(NUM_OF_PERMUTATIONS),
 		current_b2_y(NUM_OF_PERMUTATIONS), current_b3_y(NUM_OF_PERMUTATIONS),
		current_e2_y(NUM_OF_PERMUTATIONS), current_e3_y(NUM_OF_PERMUTATIONS);
	
	clock_t tt_temp=clock();
	clock_t ttt_temp;
	clock_t t_tranverse=0;
	clock_t t_cal_bounds=0;
	clock_t t_update=0;	
	
	int vertex_parent[MAX_NUM_OF_SNPS];
	
	int current_vertex,current_p,current_q,current_o,current_g,current_h,current_l;
	int count=0;
	int max_size=0;	
	
	check_points.push_back(new_root);
	//cp_a2_y.push_back(a2_y);
	//cp_a3_y.push_back(a3_y);
	
	cp_b2_y.push_back(b2_y);
	cp_b3_y.push_back(b3_y);
	
	cp_e2_y.push_back(e2_y);
	cp_e3_y.push_back(e3_y);
	
	//cp_p.push_back(p);
	cp_q.push_back(q);
	cp_o.push_back(o);
	//cp_g.push_back(g);
	cp_h.push_back(h);
	cp_l.push_back(l);
	
	vertex_parent[new_root]=-1;
	//cout <<"!!!!!";
	
	// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
	// 	cout <<a2_y[j] <<" ";
	// cout <<endl;	
	
	
	while(check_points.size()>0){
		if (max_size<check_points.size())
			max_size=check_points.size();		
        
		// count++;		
		// cout <<count <<": ";		
		// for (vector<int>::iterator iter2=check_points.begin() ; iter2!=check_points.end() ; ++iter2 ){
		// 	cout <<*iter2 <<" ";			
		// }
		// cout <<endl;
		
		current_vertex=check_points.back();		
		check_points.pop_back();
		
		// cout <<"!!!!!===";
		
		// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
		// 	cout <<current_a2_y[j] <<" ";
		// cout <<endl;
		
		//cp_a2_y.pop_back();

		// cout <<"!!!!!---";
		
		// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
		// 	cout <<current_a2_y[j] <<" ";
		// cout <<endl;


		//current_a2_y=cp_a2_y.back();
		//cp_a2_y.pop_back();
		//current_a3_y=cp_a3_y.back();
		//cp_a3_y.pop_back();
		
		current_b2_y=cp_b2_y.back();
		cp_b2_y.pop_back();
		current_b3_y=cp_b3_y.back();
		cp_b3_y.pop_back();
		
		current_e2_y=cp_e2_y.back();
		cp_e2_y.pop_back();
		current_e3_y=cp_e3_y.back();
		cp_e3_y.pop_back();		

		//current_p=cp_p.back();
		//cp_p.pop_back();

		current_q=cp_q.back();
		cp_q.pop_back();

		current_o=cp_o.back();
		cp_o.pop_back();

		//current_g=cp_g.back();
		//cp_g.pop_back();

		current_h=cp_h.back();
		cp_h.pop_back();
		
		current_l=cp_l.back();
		cp_l.pop_back();
		
		// cout <<"a2=";		
		// for (vector<int>::iterator iter=current_a2_y.begin() ; iter!=current_a2_y.end() ; ++iter){
		// 	cout <<*iter <<" ";			
		// }
		// cout <<endl;

		//int old_p=current_p;
		int old_q=current_q;
		int old_o=current_o;
		//int old_g=current_g;
		int old_h=current_h;
		int old_l=current_l;

        //cout <<"old_p=" <<old_p <<" old_g=" <<old_g <<endl;        
		//vector<int> old_a2_y(current_a2_y);
		//vector<int> old_a3_y(current_a3_y);
		vector<short> old_b2_y(current_b2_y);
		vector<short> old_b3_y(current_b3_y);
		vector<short> old_e2_y(current_e2_y);
		vector<short> old_e3_y(current_e3_y);
		
		
		for (set<edge>::iterator iter=vertex_edges[current_vertex].begin() ;
			 iter!=vertex_edges[current_vertex].end() ;
			 ++iter){			
			int other_vertex=iter->p1==current_vertex?iter->p2:iter->p1;

            //cout <<"X" <<root <<"  X" <<other_vertex <<":" <<endl;            
			//To avoid iteration, do not push the parent again!
			if (other_vertex!=vertex_parent[current_vertex]){
                
				check_points.push_back(other_vertex);				
				//current_p=old_p;
				current_q=old_q;
				current_o=old_o;
				//current_g=old_g;
				current_h=old_h;
				current_l=old_l;
				
				//current_a2_y=old_a2_y;
				//current_a3_y=old_a3_y;
				
				current_b2_y=old_b2_y;
				current_b3_y=old_b3_y;
				
				current_e2_y=old_e2_y;
				current_e3_y=old_e3_y;
				
				vertex_parent[other_vertex]=current_vertex;

				ttt_temp=clock();
				
				pair<int, int> dst_pair=make_pair(iter->p1,iter->p2);				
				SNPs_Diff & dst_diff=snps_diff[dst_pair];
				
				int sign=iter->p1==current_vertex?1:-1; //The direction of the edge will affect zero2one and one2zero.

                real_weight+=iter->weight;                
                /*
				cout <<"-->" <<"(" <<iter->p1 <<"," <<iter->p2 <<"):" <<endl;
                
                cout <<X[4].zero  <<endl;
                cout <<X[4].one  <<endl;
                cout <<X[4].two  <<endl;
                
                cout <<X[5].zero <<endl;
                cout <<X[5].one  <<endl;
                cout <<X[5].two  <<endl;

                cout <<"0->1" <<endl;                
				for (set<int>::iterator iter=dst_diff.zero2one.begin() ; iter!=dst_diff.zero2one.end() ; iter++){
                    cout <<*iter <<" ";					
				}
				cout <<endl;

                cout <<"1->0" <<endl;                
				for (set<int>::iterator iter=dst_diff.one2zero.begin() ; iter!=dst_diff.one2zero.end() ; iter++){
                    cout <<*iter <<" ";					
				}
				cout <<endl;

                cout <<"0->2" <<endl;                
				for (set<int>::iterator iter=dst_diff.zero2two.begin() ; iter!=dst_diff.zero2two.end() ; iter++){
                    cout <<*iter <<" ";					
				}
				cout <<endl;

                
                cout <<"2->0" <<endl;                
				for (set<int>::iterator iter=dst_diff.two2zero.begin() ; iter!=dst_diff.two2zero.end() ; iter++){
                    cout <<*iter <<" ";					
				}
				cout <<endl;

                cout <<"1->2" <<endl;                
				for (set<int>::iterator iter=dst_diff.one2two.begin() ; iter!=dst_diff.one2two.end() ; iter++){
                    cout <<*iter <<" ";					
				}
				cout <<endl;

                
                cout <<"2->1" <<endl;                
				for (set<int>::iterator iter=dst_diff.two2one.begin() ; iter!=dst_diff.two2one.end() ; iter++){
                    cout <<*iter <<" ";					
				}
				cout <<endl;
				*/
				//cout <<"other" <<other_vertex <<endl;

				// cout <<endl <<"!!!!!((((((";
				
				// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
				// 	cout <<current_a2_y[j] <<" ";
				// cout <<endl;

				// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
				// 	cout <<current_b2_y[j] <<" ";
				// cout <<endl;

				// cout <<"p=" <<current_p <<"," <<"q=" <<current_q <<endl;

				int * plus_data=plus_map[dst_pair];
				/*
				for (int k=0 ; k<NUM_OF_PERMUTATIONS ; ++k){
					plus_data[k]=0;					
				}
				*/
				//zero2one case				
				for (set<short>::iterator set_iter=dst_diff.zero2one.begin() ; set_iter!=dst_diff.zero2one.end() ; set_iter++){                    
					// for (vector<int>::iterator vec_iter=A_tran[*set_iter].begin() ; vec_iter!=A_tran[*set_iter].end() ; vec_iter++){
					//  	current_a2_y[*vec_iter]+=(sign);
					//  	ideal_weight++;
					//  }					
					for (vector<short>::iterator vec_iter=B_tran[*set_iter].begin() ; vec_iter!=B_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]+=(sign);						
						current_b2_y[*vec_iter]+=(sign);
						ideal_weight++;
					}
					for (vector<short>::iterator vec_iter=E_tran[*set_iter].begin() ; vec_iter!=E_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]+=(sign);						
						current_e2_y[*vec_iter]+=(sign);
						ideal_weight++;
					}
/*					
					if (X[root].zero.test(*set_iter)){
						current_p+=(sign);
					}
					else if (X[root].one.test(*set_iter)){
						current_q+=(sign);						
					}
					else{
						current_o+=(sign);						
					}
*/                    

					if (X[root].one.test(*set_iter)){
						current_q+=(sign);						
					}
					else if (X[root].two.test(*set_iter)){
						current_o+=(sign);						
					}
                    
				}//set				
				
				//one2zero case
				for (set<short>::iterator set_iter=dst_diff.one2zero.begin() ; set_iter!=dst_diff.one2zero.end() ; set_iter++){					
					// for (vector<int>::iterator vec_iter=A_tran[*set_iter].begin() ; vec_iter!=A_tran[*set_iter].end() ; vec_iter++){
					// 	current_a2_y[*vec_iter]-=(sign);
					// 	ideal_weight++;						
					// }
					for (vector<short>::iterator vec_iter=B_tran[*set_iter].begin() ; vec_iter!=B_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]-=(sign);
						current_b2_y[*vec_iter]-=(sign);
						ideal_weight++;
					}
					for (vector<short>::iterator vec_iter=E_tran[*set_iter].begin() ; vec_iter!=E_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]-=(sign);
						current_e2_y[*vec_iter]-=(sign);
						ideal_weight++;
					}
/*					
					if (X[root].zero.test(*set_iter)){
						current_p-=(sign);
					}
					else if(X[root].one.test(*set_iter)){
						current_q-=(sign);						
					}
					else {
						current_o-=(sign);						
					}
*/                    
					
					if(X[root].one.test(*set_iter)){
						current_q-=(sign);						
					}
					else if(X[root].two.test(*set_iter)){
						current_o-=(sign);						
					}

				}//set

				//zero2two case				
				for (set<short>::iterator set_iter=dst_diff.zero2two.begin() ; set_iter!=dst_diff.zero2two.end() ; set_iter++){					
					// for (vector<int>::iterator vec_iter=A_tran[*set_iter].begin() ; vec_iter!=A_tran[*set_iter].end() ; vec_iter++){
					// 	current_a3_y[*vec_iter]+=(sign);
					// 	ideal_weight++;
					// }					
					for (vector<short>::iterator vec_iter=B_tran[*set_iter].begin() ; vec_iter!=B_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]+=(sign);						
						current_b3_y[*vec_iter]+=(sign);
						ideal_weight++;
					}
					for (vector<short>::iterator vec_iter=E_tran[*set_iter].begin() ; vec_iter!=E_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]+=(sign);						
						current_e3_y[*vec_iter]+=(sign);
						ideal_weight++;
					}
/*					
					if (X[root].zero.test(*set_iter)){
						current_g+=(sign);
					}
					else if (X[root].one.test(*set_iter)){
						current_h+=(sign);						
					}
					else{
						current_l+=(sign);						
					}
*/                    

					if (X[root].one.test(*set_iter)){
						current_h+=(sign);						
					}
					else if (X[root].two.test(*set_iter)){
						current_l+=(sign);                        
					}
                    
				}//set				
				
				//two2zero case
				for (set<short>::iterator set_iter=dst_diff.two2zero.begin() ; set_iter!=dst_diff.two2zero.end() ; set_iter++){					
 					// for (vector<int>::iterator vec_iter=A_tran[*set_iter].begin() ; vec_iter!=A_tran[*set_iter].end() ; vec_iter++){
 					// 	current_a3_y[*vec_iter]-=(sign);
 					// 	ideal_weight++;						
   					// }
					for (vector<short>::iterator vec_iter=B_tran[*set_iter].begin() ; vec_iter!=B_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]-=(sign);
						current_b3_y[*vec_iter]-=(sign);
						ideal_weight++;
					}
					for (vector<short>::iterator vec_iter=E_tran[*set_iter].begin() ; vec_iter!=E_tran[*set_iter].end() ; vec_iter++){
						//plus_data[*vec_iter]-=(sign);
						current_e3_y[*vec_iter]-=(sign);
						ideal_weight++;
					}
/*					
					if (X[root].zero.test(*set_iter)){
						current_g-=(sign);
					}
					else if(X[root].one.test(*set_iter)){
						current_h-=(sign);						
					}
					else {
						current_l-=(sign);						
					}
*/                    

					if (X[root].one.test(*set_iter)){
						current_h-=(sign);						
					}
					else if (X[root].two.test(*set_iter)){
						current_l-=(sign);                        
					}
                    
				}//set

				//one2two case				
				for (set<short>::iterator set_iter=dst_diff.one2two.begin() ; set_iter!=dst_diff.one2two.end() ; set_iter++){					
  					// for (vector<int>::iterator vec_iter=A_tran[*set_iter].begin() ; vec_iter!=A_tran[*set_iter].end() ; vec_iter++){
  					// 	current_a2_y[*vec_iter]-=(sign);
  					// 	current_a3_y[*vec_iter]+=(sign);
  					// 	ideal_weight++;
  					// }					
					for (vector<short>::iterator vec_iter=B_tran[*set_iter].begin() ; vec_iter!=B_tran[*set_iter].end() ; vec_iter++){						
						//plus_data[*vec_iter]+=(sign);
						current_b2_y[*vec_iter]-=(sign);
						current_b3_y[*vec_iter]+=(sign);
						ideal_weight++;
					}
					for (vector<short>::iterator vec_iter=E_tran[*set_iter].begin() ; vec_iter!=E_tran[*set_iter].end() ; vec_iter++){						
						//plus_data[*vec_iter]+=(sign);
						current_e2_y[*vec_iter]-=(sign);
						current_e3_y[*vec_iter]+=(sign);
						ideal_weight++;
					}
/*					
					if (X[root].zero.test(*set_iter)){
						current_p-=(sign);
						current_g+=(sign);
					}
					else if (X[root].one.test(*set_iter)){
						current_q-=(sign);						
						current_h+=(sign);						
					}
					else{
						current_o-=(sign);						
						current_l+=(sign);						
					}
*/                    

					if (X[root].one.test(*set_iter)){
						current_q-=(sign);						
						current_h+=(sign);						
					}
					else if (X[root].two.test(*set_iter)){
						current_o-=(sign);						
						current_l+=(sign);                        
					}
                    
				}//set				

				//two2one case				
				for (set<short>::iterator set_iter=dst_diff.two2one.begin() ; set_iter!=dst_diff.two2one.end() ; set_iter++){					
					// for (vector<int>::iterator vec_iter=A_tran[*set_iter].begin() ; vec_iter!=A_tran[*set_iter].end() ; vec_iter++){
					// 	current_a2_y[*vec_iter]+=(sign);
					// 	current_a3_y[*vec_iter]-=(sign);
					// 	ideal_weight++;
					// }					
					for (vector<short>::iterator vec_iter=B_tran[*set_iter].begin() ; vec_iter!=B_tran[*set_iter].end() ; vec_iter++){						
						//plus_data[*vec_iter]+=(sign);
						current_b2_y[*vec_iter]+=(sign);
						current_b3_y[*vec_iter]-=(sign);
						ideal_weight++;
					}
					for (vector<short>::iterator vec_iter=E_tran[*set_iter].begin() ; vec_iter!=E_tran[*set_iter].end() ; vec_iter++){						
						//plus_data[*vec_iter]+=(sign);
						current_e2_y[*vec_iter]+=(sign);
						current_e3_y[*vec_iter]-=(sign);
						ideal_weight++;
					}
/*					
					if (X[root].zero.test(*set_iter)){
						current_p+=(sign);
						current_g-=(sign);
					}
					else if (X[root].one.test(*set_iter)){
						current_q+=(sign);						
						current_h-=(sign);						
					}
					else{
						current_o+=(sign);						
						current_l-=(sign);						
					}
*/                    

					if (X[root].one.test(*set_iter)){
						current_q+=(sign);						
						current_h-=(sign);						
					}
					else if (X[root].two.test(*set_iter)){
						current_o+=(sign);						
						current_l-=(sign);						
					}
                    
				}//set
				
				/*
				for (int k=0 ; k<NUM_OF_PERMUTATIONS ; ++k){
					current_b2_y[k]+=(plus_data[k]);					
				}

				for (int k=0 ; k<NUM_OF_PERMUTATIONS ; ++k){
					plus_data[k]=(sign*plus_data[k]);					
				}
				*/
				
  				
				for (int kk=0 ; kk<NUM_OF_PERMUTATIONS ; ++kk){
                    //cout << other_vertex <<":";                    
                    // cout <<XiY_B[other_vertex][kk] <<"="
                    //      <<current_a2_y[kk] <<" + "
                    //      <<current_b2_y[kk] <<" + "
                    //      <<current_e2_y[kk] <<endl;
                    
					current_a2_y[kk]=XiY_B[other_vertex][kk]-current_b2_y[kk]-current_e2_y[kk];
					current_a3_y[kk]=XiY_E[other_vertex][kk]-current_b3_y[kk]-current_e3_y[kk];
				}
                
				
				current_p=X[other_vertex].one.count()-current_q-current_o;
				current_g=X[other_vertex].two.count()-current_h-current_l;
				
					
				t_update+=(clock()-ttt_temp);				
				//cout <<"!!!!!]]]]]]]]";
				
				// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
				// 	cout <<current_a2_y[j] <<" ";
				// cout <<endl;

				// for (int j=0; j<NUM_OF_PERMUTATIONS ; ++j)
				// 	cout <<current_b2_y[j] <<" ";
				// cout <<endl;
				
				//cout <<"af p=" <<current_p <<" " <<"g=" <<current_g <<endl;
				//cout <<diff.zero2one.size() <<" " <<diff.one2zero.size() <<endl;


				ttt_temp=clock();				
				calculate_test_scores(data_map,
                                      current_a2_y,current_a3_y,
                                      current_b2_y,current_b3_y,
                                      current_e2_y,current_e3_y,
                                      X[root].zero.count()-current_p-current_g,
                                      X[root].one.count()-current_q-current_h,
                                      X[root].two.count()-current_o-current_l,
                                      current_p,current_q,current_o,
                                      current_g,current_h,current_l,
                                      abcdef,bins,org_bins,chi_alpha,chi_alpha_array);
				
				t_cal_bounds+=(clock()-ttt_temp);				
					
				//cp_a2_y.push_back(current_a2_y);
				//cp_a3_y.push_back(current_a3_y);
				
				cp_b2_y.push_back(current_b2_y);
				cp_b3_y.push_back(current_b3_y);
				
				cp_e2_y.push_back(current_e2_y);
				cp_e3_y.push_back(current_e3_y);
				
				//cp_p.push_back(current_p);
				cp_q.push_back(current_q);
				cp_o.push_back(current_o);
				//cp_g.push_back(current_g);
				cp_h.push_back(current_h);
				cp_l.push_back(current_l);
			}
		}
	}
	
	// cout <<"Max check points size:" <<max_size <<endl;
	tt_traverse_tree+=(clock()-tt_temp-t_cal_bounds-t_update);
	tt_update_old+=(t_update);
	tt_cal_bounds+=(t_cal_bounds);
	old_update_times++;	
}



void * roam_tree(void *thread_id)
{
    int *abcdef[6]={0};
    double chi_alpha=0;

    double *chi_alpha_array;
    
    __int64 visited_points=0;
    __int64 total_points=0;

    __int64 *bins;
    __int64 *org_bins;    
    
    __int64 ideal_weight=0;
    __int64 real_weight=0;
    
    clock_t tt_init=0,tt_cal_bounds=0,tt_visit_tree=0,tt_traverse_tree=0,tt_cut_tree=0,tt_update=0;
    clock_t tt_update_new=0;
    clock_t tt_update_old=0;
    clock_t tt_write_file=0;

    __int64 update_rows=0;
    __int64 new_update_times=0;
    __int64 old_update_times=0;

    
    clock_t t_tmp;	
    //* Initialize Bins *//
    bins=new __int64[n_bins];
    memset(bins, 0, n_bins*sizeof(__int64));

    org_bins=new __int64[n_bins];
    memset(org_bins, 0, n_bins*sizeof(__int64));


    chi_alpha_array=new double[NUM_OF_PERMUTATIONS];
    memset(chi_alpha_array, 0 ,NUM_OF_PERMUTATIONS*sizeof(double));    
        
    long id=(long)thread_id;    
    //*Roaming on MST*//	
	clock_t tt_temp;
	
	cout <<"[" <<id <<"] " <<"Roaming on MST..." <<endl;
	cout.flush();	
	t_tmp = clock();
	//long real_mst_weight=0;

	//Cut the max weight one degree edge
	//set<edge, greater<edge> > leave_edge;
	
	//Cut the min weight one degree edge
	set<edge > leave_edge;
	
	vector< set<edge> > vertex_edges(MAX_NUM_OF_SNPS);

	tt_temp=clock();	
	for (list<edge>::iterator iter=mst.begin() ; iter!=mst.end() ; iter++){
		vertex_edges[iter->p1].insert(*iter);
		vertex_edges[iter->p2].insert(*iter);
	}
	
	// for (set<edge>::iterator iter=vertex_edges[21].begin() ; iter!=vertex_edges[21].end() ; ++iter){
	// 	cout <<"(" <<iter->p1 <<"," <<iter->p2 <<"), ";
	// }
	// cout <<endl;	
	
	for (int jj=0; jj<MAX_NUM_OF_SNPS ; ++jj){
		if (vertex_edges[jj].size()==1){
			// cout <<jj <<") ";			
			for (set<edge>::iterator iter=vertex_edges[jj].begin() ; iter!=vertex_edges[jj].end() ; ++iter){
				// cout <<"(" <<iter->p1 <<"," <<iter->p2 <<"), ";
				leave_edge.insert(*iter);
			}
			// cout <<endl;			
		}		
	}	
	// cout <<endl;

	#ifdef verbose
	cout <<"Initial one degree edges:";	
	for (set<edge>::iterator iter2=leave_edge.begin() ; iter2!=leave_edge.end() ;++iter2){
	 	cout <<"(" <<iter2->p1 <<"," <<iter2->p2 <<"), ";
	 }
	cout <<endl;
	#endif
	
	tt_cut_tree+=(clock()-tt_temp);

    abcdef[0]=new int[NUM_OF_PERMUTATIONS];
	abcdef[1]=new int[NUM_OF_PERMUTATIONS];
	abcdef[2]=new int[NUM_OF_PERMUTATIONS];
	abcdef[3]=new int[NUM_OF_PERMUTATIONS];
	abcdef[4]=new int[NUM_OF_PERMUTATIONS];
	abcdef[5]=new int[NUM_OF_PERMUTATIONS];
    
	
	map< pair<int,int>, int* > previous_plus;
/*    
	for (list<edge>::iterator iter=mst.begin() ; iter!=mst.end() ; ++iter){
		pair<int,int> point_pair=make_pair(iter->p1,iter->p2);
		int * plus_data;
		plus_data=new int[NUM_OF_PERMUTATIONS];
		for (int i=0 ; i<NUM_OF_PERMUTATIONS ; ++i)
			plus_data[i]=0;		
		previous_plus.insert(make_pair(point_pair,plus_data));		
	}
*/
	vector< vector<int> > previous_B_transpose(MAX_NUM_OF_INDIVIDUALS);

	int previous_node=-1;
	
	int next_cut_node=-1;	
	set<edge>::iterator next_cut_leaf=leave_edge.end();


    __int64 from,to,size;
    
    if (!get_snps_partition(MAX_NUM_OF_SNPS,n_partitions,id,from,to,size)) {
        cout <<"No work!" <<endl;
        //exit(1);        
        //return 1;
        pthread_exit(NULL);        
    }
    
    __int64 snps_range_size=to-from+1;
    __int64 o_job_percent=0;
    
    cout <<"[" <<id <<"] " <<"Parallel SNP range:" <<from <<"--" <<to <<"\t: " <<size <<" pairs" <<endl;
	for (int jj=0; jj<(MAX_NUM_OF_SNPS-1) ; ++jj){		
		set<edge>::iterator temp;
		tt_temp=clock();

		if (next_cut_leaf==leave_edge.end()){			
			if (leave_edge.size()>0){
				temp=leave_edge.begin();
			}
			else {
				cout <<"Unknown error!" <<endl;			
				break;
			}
		}
		else {
			temp=next_cut_leaf;
		}
		int p1=temp->p1;
		int p2=temp->p2;
		int root=-1;
		int check_vertex=-1;		

		// cout <<p1 <<"--";
		
		// for (set<edge>::iterator iter=vertex_edges[p1].begin() ; iter!=vertex_edges[p1].end() ; ++iter){
		// 	cout <<iter->p1 <<" " <<iter->p2 <<",";
			
		// }
		// cout <<endl;

		// cout <<p2 <<"--";
		
		// for (set<edge>::iterator iter=vertex_edges[p2].begin() ; iter!=vertex_edges[p2].end() ; ++iter){
		// 	cout <<iter->p1 <<" " <<iter->p2 <<",";
			
		// }
		// cout <<endl;		
		if (vertex_edges[p1].size()==1) {
			root=p1;
			check_vertex=p2;			
		}
		else {
			root=p2;
			check_vertex=p1;			
		}
		
		vertex_edges[p1].erase(*temp);
		vertex_edges[p2].erase(*temp);

		leave_edge.erase(temp);
		//cout <<endl <<"--Cut (" <<root <<"," <<check_vertex <<")" <<endl;
		if (vertex_edges[check_vertex].size()==1){			
/*			for (set<edge>::iterator iter=vertex_edges[check_vertex].begin() ; iter!=vertex_edges[check_vertex].end() ; ++iter){
				leave_edge.insert(*iter);
			}
*/
			set<edge>::iterator iter=vertex_edges[check_vertex].begin();
			pair<set<edge>::iterator,bool> ret;
			ret=leave_edge.insert(*(iter));
//			if (ret.second==false) {				
//				cout <<"Unknown error!" <<endl;
//				break;
//			}
//			else{
				next_cut_leaf=ret.first;
				next_cut_node=check_vertex;				
//			}		
            #ifdef verbose				
			cout <<"---Set next leaf to be cut: ("
				 <<next_cut_leaf->p1 <<"," <<next_cut_leaf->p2 <<")"
				 <<endl;
			#endif
		}
		else {
			next_cut_leaf=leave_edge.end();
			next_cut_node=-1;			
		}		
		// cout <<jj <<" " <<"(" <<temp->p1 <<"," <<temp->p2 <<"):" <<temp->weight <<endl;
		
		#ifdef verbose
		cout <<endl <<"Cut (" <<root <<"," <<check_vertex <<")" <<endl;
		cout <<"One degree edges after cutting:";
		for (set<edge>::iterator iter2=leave_edge.begin() ; iter2!=leave_edge.end() ;++iter2){
			cout <<"(" <<iter2->p1 <<"," <<iter2->p2 <<"), ";
		}
		cout <<endl;

		#endif
			
		tt_cut_tree+=(clock()-tt_temp);

        if (jj<from)
            continue;
        else if (jj>to)
            break;

        real_weight+=(temp->weight);

		ideal_weight+=temp->weight*Ylist.size();		

/*Initialize updating seed*/

		tt_temp=clock();
		triSNP x1,x2;

		x1=X[root];
		x2=X[check_vertex];
		
		vector<Intermedia> A2_Y(NUM_OF_PERMUTATIONS);
		vector<short> a2_y(NUM_OF_PERMUTATIONS);
		vector<Intermedia> A3_Y(NUM_OF_PERMUTATIONS);
		vector<short> a3_y(NUM_OF_PERMUTATIONS);

		vector<Intermedia> B2_Y(NUM_OF_PERMUTATIONS);
		vector<short> b2_y(NUM_OF_PERMUTATIONS);	
		vector<Intermedia> B3_Y(NUM_OF_PERMUTATIONS);
		vector<short> b3_y(NUM_OF_PERMUTATIONS);
		
		vector<Intermedia> E2_Y(NUM_OF_PERMUTATIONS);
		vector<short> e2_y(NUM_OF_PERMUTATIONS);
		vector<Intermedia> E3_Y(NUM_OF_PERMUTATIONS);
		vector<short> e3_y(NUM_OF_PERMUTATIONS);
		
		
		biPhenotype y;
		int a,b,c,d,e,f,p,q,o,g,h,l,s,t,r;
		Intermedia C,D,F,P,Q,O,G,H,L;
		
		vector<Intermedia> A(NUM_OF_PERMUTATIONS),B(NUM_OF_PERMUTATIONS),E(NUM_OF_PERMUTATIONS);
		
		Data_map data_map;

		vector< vector<short> > A_transpose(MAX_NUM_OF_INDIVIDUALS), B_transpose(MAX_NUM_OF_INDIVIDUALS), E_transpose(MAX_NUM_OF_INDIVIDUALS);		
			
		list<biPhenotype>::iterator yp;

        //cout <<"Calculate P" <<endl;        
		P=x1.zero & x2.one;
		p=P.count();			
		G=x1.zero & x2.two;
		g=G.count();
		s=x1.zero.count()-p-g;		

		Q=x1.one & x2.one;
		q=Q.count();			
		H=x1.one & x2.two;
		h=H.count();
		t=x1.one.count()-q-h;		

		O=x1.two & x2.one;
		o=O.count();
		L=x1.two & x2.two;
		l=L.count();
		r=x1.two.count()-o-l;		

		//Build the entry map
		int yp_count=0;		
		for (yp=Ylist.begin() ; yp!=Ylist.end() ; yp++) {
			/*Calculate A,B,C,D,E,F*/				
			A[yp_count]=x1.zero & ~(*yp);
			a=A[yp_count].count();
			abcdef[0][yp_count]=a;			
				
			B[yp_count]=x1.one & ~(*yp);
			b=B[yp_count].count();
			abcdef[1][yp_count]=b;
			
			E[yp_count]=x1.two & ~(*yp);
			e=E[yp_count].count();
			abcdef[4][yp_count]=e;
			
			C=x1.zero & (*yp);				
			c=C.count();
			abcdef[2][yp_count]=c;
			  
			D=x1.one & (*yp);
			d=D.count();
			abcdef[3][yp_count]=d;			

			F=x1.two & (*yp);
			f=F.count();
			abcdef[5][yp_count]=f;
			
			for (int o=0 ; o<MAX_NUM_OF_INDIVIDUALS ; ++o){
				
				if (A[yp_count].test(o)){
					A_transpose[o].push_back(yp_count);					
				}				
				else if(B[yp_count].test(o)){
					B_transpose[o].push_back(yp_count);					
				}
				else if(E[yp_count].test(o)){
					E_transpose[o].push_back(yp_count);					
				}
				/*
				if(B[yp_count].test(o)){
					B_transpose[o].push_back(yp_count);					
				}
				*/
			}
			
			//cout <<endl <<A <<endl <<a <<endl <<B <<endl <<b <<endl ;
			//cout <<endl <<P <<endl <<p <<endl <<Q <<endl <<q <<endl ;
			
/*			
			pair <int,int> index;
			index=make_pair(b,e);	
			//cout <<index.first <<" " <<index.second <<endl;
		
			Data_map::iterator pos;
			pos = data_map.find(index);

			if (pos == data_map.end()){			
				vector<int> tmp_vector;
				tmp_vector.push_back(yp_count);			
				data_map.insert(make_pair(index,tmp_vector));
			}
			else {
				data_map[index].push_back(yp_count);			
			}
*/			
			
			A2_Y[yp_count]=A[yp_count] & P;
			A3_Y[yp_count]=A[yp_count] & G;

			B2_Y[yp_count]=B[yp_count] & Q;
			B3_Y[yp_count]=B[yp_count] & H;

			E2_Y[yp_count]=E[yp_count] & O;
			E3_Y[yp_count]=E[yp_count] & L;			
		
			// old_range[yp_count][0]=max(0,p-c);
			// old_range[yp_count][1]=min(p,a);
		
			// old_range[yp_count][2]=max(0,q-d);
			// old_range[yp_count][3]=min(q,b);
          
			// cout <<"p=" <<p <<"a=" <<a <<"c=" <<c
			// <<"q=" <<q <<"b=" <<b <<"d=" <<d
			// <<endl;

			//cerr <<a <<"\t" <<b <<endl;			
			yp_count++;
		}

		//cout <<"Total Num of Entries:" <<data_map.size() <<endl;
	
		for (int l=0 ; l<NUM_OF_PERMUTATIONS ; ++l){
			a2_y[l]=A2_Y[l].count();
			a3_y[l]=A3_Y[l].count();
			
			b2_y[l]=B2_Y[l].count();
			b3_y[l]=B3_Y[l].count();
			
			e2_y[l]=E2_Y[l].count();
			e3_y[l]=E3_Y[l].count();
		}
		tt_init+=(clock()-tt_temp);
		
		//cout <<"Initial point" <<endl;
		tt_temp=clock();		
		calculate_test_scores(data_map,a2_y,a3_y,b2_y,b3_y,e2_y,e3_y,s,t,r,p,q,o,g,h,l,abcdef,bins,org_bins,chi_alpha,chi_alpha_array);
		tt_cal_bounds+=(clock()-tt_temp);		
		//cout <<endl;
		#ifdef verbose
		for (Data_map::iterator iter=data_map.begin() ; iter!=data_map.end() ; iter++){
			cout <<iter->first.first <<" " <<iter->first.second <<":";
			for(vector<int>::iterator iter2=iter->second.begin() ; iter2!=iter->second.end() ; iter2++){
				cout <<(*iter2) <<" ";			
			}
			cout <<endl;		
		}
		#endif
		
/////////////////////	
/*Update*/
	
		//cout <<endl <<"Traverse begin" <<endl;
		//traverse_mst(root,check_vertex,vertex_edges,snps_diff,data_map,X,XiY_B,a2_y,b2_y,A_transpose,B_transpose,previous_B_transpose,previous_plus,previous_node,p,q);
		
//		if (previous_node==-1)
			//A new one degree node in a different lineage, recalculate B_transpose
		
			traverse_mst(root,check_vertex,vertex_edges,snps_diff,data_map,X,
						 XiY_B,XiY_E,
						 a2_y,a3_y,b2_y,b3_y,e2_y,e3_y,
						 A_transpose,B_transpose,E_transpose,
						 previous_plus,
						 p,q,o,g,h,l,abcdef,chi_alpha,chi_alpha_array,
                         real_weight,ideal_weight,
                         tt_traverse_tree,tt_update_old,tt_cal_bounds,old_update_times,bins,org_bins);
						 
//		else
			//A one degree node in the same lineage, use previous_B_transpose to improve performance
//			traverse_mst(root,check_vertex,vertex_edges,snps_diff,data_map,X,XiY_B,a2_y,b2_y,A_transpose,B_transpose,previous_B_transpose,previous_plus,previous_node,p,q);
		//cout <<"Traverse end" <<endl <<endl;
/*
		if (next_cut_node!=-1){			
			previous_node=root;
			previous_B_transpose=B_transpose;
		}
		else {
			previous_node=-1;
		}
*/
/*
		if ((jj+1)%100==0) {
            cout <<"[" <<id <<"] " <<jj+1 <<"/" <<(to) <<": " <<chi_alpha <<endl;            
		}
*/
            
            __int64 job_percent=(((jj-from+1)*100/snps_range_size));
            
            if ((job_percent-o_job_percent)>=(1)){
                o_job_percent=job_percent;
                cout <<"[" <<id <<"] " <<setw(5) <<job_percent/1 <<"% : " <<chi_alpha <<endl;
                cout <<"@@" <<int(job_percent*0.9) <<endl;                
            }            
	}
/*    
    if ((jj+1)%100) {
        cout <<"[" <<id <<"] " <<jj+1 <<"/" <<(to) <<": " <<chi_alpha <<endl;            
    }
*/    
	cout <<endl  <<"[" <<id <<"] " <<setw(5) <<"100" <<"% : " <<"Max Test Score=" <<chi_alpha <<endl <<endl;

    cout <<"[" <<id <<"] " <<"Max Test Score original phenotypes:"
         <<chi_alpha_array[0]
         <<endl <<endl;
    
    cout <<"[" <<id <<"] " <<"Max Test Score for each permutations:";    
    for (int i=1 ; i<NUM_OF_PERMUTATIONS ; ++i)        
        cout <<chi_alpha_array[i] <<"\t";    
    cout <<endl <<endl;
    
	cout <<"[" <<id <<"] " <<"Real mst weight:" <<real_weight <<endl;
	cout <<"[" <<id <<"] " <<"(B) Percentage:" <<(real_weight*100.0/size)/MAX_NUM_OF_INDIVIDUALS <<"%" <<endl;
	
	cout <<"[" <<id <<"] " <<"Time Elapsed for Visit tree(Cut & traverse):" <<(tt_cut_tree+tt_traverse_tree)*1.0/CLOCKS_PER_SEC <<" s"
        //<< tt_cut_tree <<"\t" <<tt_traverse_tree
         <<endl;
	cout <<"[" <<id <<"] " <<"Time Elapsed for Calculate Test Score:" <<(tt_cal_bounds)*1.0/CLOCKS_PER_SEC <<" s" <<endl;
	cout <<"[" <<id <<"] " <<"Time Elapsed for Update:" <<(tt_init+tt_update_old+tt_update_new)*1.0/CLOCKS_PER_SEC <<" s" <<endl;
	//cout <<"[" <<endl;
	//cout <<"[" <<id <<"] " <<"Time Elapsed for Init:" <<(tt_init)*1.0/CLOCKS_PER_SEC <<" s" <<endl;
	//cout <<"[" <<id <<"] " <<"Time Elapsed for Update_old:" <<(tt_update_old)*1.0/CLOCKS_PER_SEC <<" s" <<endl;
	//cout <<"[" <<id <<"] " <<"Time Elapsed for Update_new:" <<(tt_update_new)*1.0/CLOCKS_PER_SEC <<" s" <<endl;
	//cout <<"]" <<endl;
	
	cout <<endl;	

	//cout <<"[" <<id <<"] " <<"Old update run times:" <<old_update_times <<endl;
	//cout <<"[" <<id <<"] " <<"New update run times:" <<new_update_times <<endl;	
	//cout <<endl;
	
	cout <<"[" <<id <<"] " <<"Ideal mst weight:" <<ideal_weight <<endl;
	//cout <<"[" <<id <<"] " <<"Update rows:" <<update_rows <<endl;	
	cout <<"[" <<id <<"] " <<"(C) Percentage:" <<(ideal_weight*100.0/size)/NUM_OF_PERMUTATIONS/(MAX_NUM_OF_INDIVIDUALS) <<"%" <<endl;
	cout <<"[" <<id <<"] " <<"(D) Pruning ratio:" <<100-(ideal_weight*100.0/size)/NUM_OF_PERMUTATIONS/(MAX_NUM_OF_INDIVIDUALS) <<"%" <<endl;
	cout <<endl;
    
/*test*/
	/*
	Intermedia tmp;
	tmp=(~X[6].genotype) & (X[8].genotype);
	cout <<tmp.count() <<endl;
	tmp=X[6].genotype & X[8].genotype;
	cout <<tmp.count() <<endl;
	*/

    cout <<"[" <<id <<"] " <<"All Done!" <<endl <<endl;

    ofstream pvalue_output;
    pvalue_output.open(pvalue_fn);    
    pvalue(pvalue_output,bins,org_bins);
    
/*    
    ofstream bins_output;
    if (bins_fn==NULL){        
        cout.precision(10);
        //cout <<bin_interval <<endl;    
        for (int i=0 ; i<n_bins ; ++i){
            cout <<id <<"\t" <<i*bin_interval <<"\t" <<bins[i] <<"\t" <<org_bins[i] <<endl;
        }        
    }
    else {        
        char fn_suffix[10];
        sprintf(fn_suffix,".%03ld",id);
        string bins_fn_with_id(bins_fn);
        bins_fn_with_id.append(fn_suffix);
        
        bins_output.open(bins_fn_with_id.c_str());
        for (int i=0 ; i<n_bins ; ++i){
            bins_output <<i*bin_interval <<"\t" <<bins[i] <<"\t" <<org_bins[i] <<endl;
        }
        bins_output.close();        
    }
*/    
    delete []bins;
    delete []org_bins;
    delete []chi_alpha_array;    
    pthread_exit(NULL);
}

//////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	clock_t t_start = clock();
	clock_t t_tmp;	
	int i;
    
	get_args(argc,argv);

    bin_interval=MAX_NUM_OF_INDIVIDUALS*1.0/n_bins;
    
    
    //*Read Data*//
    if (read_tri_genotypes(geno_infile,X)<0){
        exit(1);        
    }    

    biPhenotype Y_origin(MAX_NUM_OF_INDIVIDUALS);
    NUM_OF_ONES_IN_Y=read_bi_phenotypes(pheno_infile,Y_origin);
    
    Ylist.push_back(Y_origin);

    //* Permutations *//
    t_tmp = clock();
	cout <<"Permutating phenotypes...";
	cout.flush();	
	generate_permutation(NUM_OF_PERMUTATIONS,NUM_OF_ONES_IN_Y,Ylist);
	cout <<"OK" <<endl;
    cout <<NUM_OF_PERMUTATIONS <<" permutations generated." <<endl;	
	clock_t t_permutation = clock();	
	cout <<"Time Elapsed for permutations:" <<(t_permutation-t_tmp)*1.0/CLOCKS_PER_SEC <<" s" <<endl <<endl;    

    ++NUM_OF_PERMUTATIONS; // Consider the original phenotypes a permutation
    
    //*Building Minimum Spanning Tree*//	
	t_tmp = clock();	
	cout <<"Building Minimum Spanning Tree...";
	cout.flush();	
	int mst_total_weight=build_mst(X,mst);
	//int mst_total_weight=build_linear_tree(X,mst);
	cout <<"OK" <<endl;	
	clock_t t_gen_mst = clock();
	
	// __int64 total_weight=0;	
	// cout <<"Total Weight for MST:" <<mst_total_weight <<endl;
	// for (int ii=0 ; ii<MAX_NUM_OF_SNPS ; ii++){
	// 	for (int jj=ii+1 ; jj< MAX_NUM_OF_SNPS ; jj++){
	// 		total_weight+=dist_matrix[ii][jj];			
	// 	}		
	// }

	#ifdef verbose
		cout <<"MST edges:";	
		print_mst(mst);
	#endif
	
	cout <<"MST size:" <<mst.size() <<endl;	
	cout <<"MST weight:" <<mst_total_weight <<endl;	
	cout <<"(A) Percentage:" <<mst_total_weight*100.0/(MAX_NUM_OF_SNPS*MAX_NUM_OF_INDIVIDUALS) <<"%" <<endl;	
	cout <<"Time Elapsed for Building MST:" <<(t_gen_mst-t_tmp)*1.0/CLOCKS_PER_SEC <<" s" <<endl <<endl;

    
    //*Calculate snps difference in mst*//	
	t_tmp = clock();	
	cout <<"Caculating SNPs Pairwise difference in MST..." <<endl;
	cout.flush();
    //cout <<a1 <<" " <<a2 <<" " <<a3 <<" " <<c1 <<" " <<c2 <<" " <<c3 <<" " <<S <<" " <<P <<" " <<G <<result <<endl;
    calculate_snps_diff(X,mst,snps_diff);
	cout <<"OK" <<endl;	
	clock_t t_cal_diff = clock();	
	cout <<"Time Elapsed for Calculate Diff:" <<(t_cal_diff-t_tmp)*1.0/CLOCKS_PER_SEC <<" s" <<endl <<endl;

    
    //*Calculate XiY_B, XiY_E*//
	t_tmp = clock();
	cout <<"Caculating XiY_B,XiY_E..." <<endl;
	cout.flush();

    XiY_B=new short*[MAX_NUM_OF_SNPS];
    
	for (int ll=0 ; ll<MAX_NUM_OF_SNPS ; ++ll){
		XiY_B[ll]=new short[NUM_OF_PERMUTATIONS];
	}
    
    XiY_E=new short*[MAX_NUM_OF_SNPS];
    
	for (int ll=0 ; ll<MAX_NUM_OF_SNPS ; ++ll){
		XiY_E[ll]=new short[NUM_OF_PERMUTATIONS];
	}


	Intermedia tmp;
	int ii,jj;

    __int64 progress,o_progress=0;
	ii=0;
	for (vector<triSNP>::iterator iter1=X.begin(); iter1!=X.end() ; ++iter1, ++ii){
		jj=0;		
		for (list <biPhenotype>::iterator iter2=Ylist.begin(); iter2!=Ylist.end() ; ++iter2, ++jj){
			tmp=(iter1->one) & (~(*iter2)); //Xi=1 Y=0
			XiY_B[ii][jj]=tmp.count();
			tmp=(iter1->two) & (~(*iter2)); //Xi=2 Y=0
			XiY_E[ii][jj]=tmp.count();
		}

        progress=(ii+1)*1000/X.size();
        if ((progress-o_progress)>=1) {
            o_progress=progress;
            
            cout <<setw(5) <<progress/10.0 <<"%";
            cout <<"\r";
            cout.flush();
        }        
	}
    cout <<endl;	
	cout <<"OK" <<endl;
	clock_t t_cal_xiy=clock();	
	cout <<"Time Elapsed for Calculate XiY_B,XiY_E:" <<(t_cal_xiy-t_tmp)*1.0/CLOCKS_PER_SEC <<" s" <<endl <<endl;
    
    
    //*Look up tables for acceleration of computing*//
    look_up_n0=new double* [MAX_NUM_OF_INDIVIDUALS+1];
    for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS+1 ; ++i){
        look_up_n0[i]=new double [MAX_NUM_OF_INDIVIDUALS+1];
    }
    look_up_n1=new double* [MAX_NUM_OF_INDIVIDUALS+1];
    for (int i=0 ; i<MAX_NUM_OF_INDIVIDUALS+1 ; ++i){
        look_up_n1[i]=new double [MAX_NUM_OF_INDIVIDUALS+1];
    }
    create_chi_square_lookup_tables(look_up_n0,look_up_n1);
    //cout <<look_up_n0[0][1];    


    pthread_t threads[n_threads];
    int rc;
    long t;
    for(t=0; t<n_threads; t++){
        printf("In main: creating thread %ld\n", id+t);
        rc = pthread_create(&threads[t], NULL, roam_tree, (void *)(id+t));
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    //Wait until the tread.
    for(t=0; t<n_threads; t++){
        pthread_join(threads[t],NULL);        
    }


    cout <<"Total Runtime:" <<(clock()-t_start)*1.0/CLOCKS_PER_SEC <<" s" <<endl <<endl;    
    
/*Recycle*/    
    for (int si=0 ; si<MAX_NUM_OF_INDIVIDUALS+1 ; ++si){
        delete []look_up_n0[si];
        delete []look_up_n1[si];
    }

    for (int ll=0 ; ll<MAX_NUM_OF_SNPS ; ++ll){
		delete []XiY_B[ll];
        delete []XiY_E[ll];
	}

    pthread_exit(NULL);
}
