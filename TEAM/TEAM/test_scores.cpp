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

#include<math.h>
#include<unistd.h>
#include<sys/types.h>
#include<stdio.h>
#include<stdlib.h>
#include<limits.h>

#include "global_param.hpp"

extern int NUM_OF_ONES_IN_Y;
extern int MAX_NUM_OF_INDIVIDUALS;
extern double **look_up_n0;
extern double **look_up_n1;

//Basic chi_square calculation. Abandoned in the current version.
double chi_square_basic(int a2, int a3, int A, int C, int P, int G, int n1){	
	int a1,c1,c2,c3,n0,S,M;
	double Ea1,Ea2,Ea3,Ec1,Ec2,Ec3;
	M=MAX_NUM_OF_INDIVIDUALS;
	

	a1=A-a2-a3;
	c2=P-a2;
	c3=G-a3;	
	c1=C-c2-c3;
	S=a1+c1;
	
	n0=M-n1;
	
	Ea1=(n0*S*1.0)/(M);
	Ea2=(n0*P*1.0)/(M);
	Ea3=(n0*G*1.0)/(M);
	
	Ec1=(n1*S*1.0)/(M);
	Ec2=(n1*P*1.0)/(M);
	Ec3=(n1*G*1.0)/(M);
	
	double result=0;
	
	if (Ea1!=0)
		result+=(a1-Ea1)*(a1-Ea1)/Ea1;
	
	if (Ea2!=0)
		result+=(a2-Ea2)*(a2-Ea2)/Ea2;

	if (Ea3!=0)
		result+=(a3-Ea3)*(a3-Ea3)/Ea3;

	if (Ec1!=0)
		result+=(c1-Ec1)*(c1-Ec1)/Ec1;
	
	if (Ec2!=0)
		result+=(c2-Ec2)*(c2-Ec2)/Ec2;

	if (Ec3!=0)
		result+=(c3-Ec3)*(c3-Ec3)/Ec3;

	//result=(a1-Ea1)*(a1-Ea1)/Ea1+(a2-Ea2)*(a2-Ea2)/Ea2+(c1-Ec1)*(c1-Ec1)/Ec1+(c2-Ec2)*(c2-Ec2)/Ec2;	
	// cout <<" A=" <<A <<" C=" <<C <<" P="  <<P <<" M=" <<M <<endl;
	// cout <<" a1=" <<a1 <<" a2=" <<a2 <<" c1=" <<c1 <<" c2="  <<c2 <<" result=" <<result <<endl <<endl;

	return result;	
}

//An improved version of chi square test computing. Abandoned in the current version as well.
inline double chi_square_improved(int a2, int a3, int A, int C, int P, int G,					   
					   double Ea1, double Ea2, double Ea3, double Ec1, double Ec2, double Ec3){
	int a1,c1,c2,c3;	
	a1=A-a2-a3;
	c2=P-a2;
	c3=G-a3;	
	c1=C-c2-c3;
	
	
	double result=0;
	
	if (Ea1!=0)
		result+=(a1-Ea1)*(a1-Ea1)/Ea1;
	
	if (Ea2!=0)
		result+=(a2-Ea2)*(a2-Ea2)/Ea2;

	if (Ea3!=0)
		result+=(a3-Ea3)*(a3-Ea3)/Ea3;

	if (Ec1!=0)
		result+=(c1-Ec1)*(c1-Ec1)/Ec1;
	
	if (Ec2!=0)
		result+=(c2-Ec2)*(c2-Ec2)/Ec2;

	if (Ec3!=0)
		result+=(c3-Ec3)*(c3-Ec3)/Ec3;

	//result=(a1-Ea1)*(a1-Ea1)/Ea1+(a2-Ea2)*(a2-Ea2)/Ea2+(c1-Ec1)*(c1-Ec1)/Ec1+(c2-Ec2)*(c2-Ec2)/Ec2;
	
	// cout <<" A=" <<A <<" C=" <<C <<" P="  <<P <<" M=" <<M <<endl;
	// cout <<" a1=" <<a1 <<" a2=" <<a2 <<" c1=" <<c1 <<" c2="  <<c2 <<" result=" <<result <<endl <<endl;

	return result;	
}

//Use a look-up table to accelerate the computing efficiency. IN USE NOW.
double chi_square(int a2, int a3, int A, int C, int S, int P, int G){    
	int a1,c1,c2,c3;	
	a1=A-a2-a3;
	c2=P-a2;
	c3=G-a3;	
	c1=C-c2-c3;
    
	double result=0;
    //cout <<"chi begin" <<endl;
    //cout <<"    " <<a1 <<" " <<a2 <<" " <<a3 <<" " <<c1 <<" " <<c2 <<" " <<c3 <<" " <<S <<" " <<P <<" " <<G <<" " <<endl;    
    result=look_up_n0[a1][S]+look_up_n0[a2][P]+look_up_n0[a3][G]+look_up_n1[c1][S]+look_up_n1[c2][P]+look_up_n1[c3][G];
    //cout <<"chi end" <<endl;
    
    //cout <<"    " <<a1 <<" " <<a2 <<" " <<a3 <<" " <<c1 <<" " <<c2 <<" " <<c3 <<" " <<S <<" " <<P <<" " <<G <" " <<result <<endl;    
    //cout <<result <<look_up_n0[a1][S] <<endl;    
    //cout <<"\t"<<result <<endl;
    //cout <<look_up_n0[a1][S] <<" " <<look_up_n0[a2][P] <<" " <<look_up_n0[a3][G] <<endl;
    //cout <<look_up_n1[c1][S] <<" " <<look_up_n1[c2][P] <<" " <<look_up_n1[c3][G] <<endl;
	return result;	
}


int create_chi_square_lookup_tables(double **look_up_n0, double **look_up_n1)
{
    cout <<"Generating Lookup Tables...";
    
    int n1=NUM_OF_ONES_IN_Y;
    int n0=MAX_NUM_OF_INDIVIDUALS-n1;    
    int M=MAX_NUM_OF_INDIVIDUALS;
    double look_up_tmp;
    int n0s,n1s;    
    for (int si=0 ; si<MAX_NUM_OF_INDIVIDUALS+1 ; ++si){
        //n0s=n0*si;
        for (int ai=0 ; ai<MAX_NUM_OF_INDIVIDUALS+1 ; ++ai){
            if (si==0){
                look_up_n0[ai][si]=0;                
            }
            else {                
                //ook_up_tmp=M*ai-n0s;            
                //look_up_n0[ai][si]=look_up_tmp*look_up_tmp/n0s;
                look_up_n0[ai][si]=(ai-n0*si*1.0/M)*(ai-n0*si*1.0/M)/(n0*si*1.0/M);                
            }            
        }        
    }

    for (int si=0 ; si<MAX_NUM_OF_INDIVIDUALS+1 ; ++si){
        //n1s=n1*si;
        //look_up_n1[0][n1s]=0;        
        for (int ai=0 ; ai<MAX_NUM_OF_INDIVIDUALS+1 ; ++ai){
            if (si==0)
                look_up_n1[ai][si]=0;
            else {                
                //look_up_tmp=M*ai-n1s;            
                //look_up_n1[ai][si]=look_up_tmp*look_up_tmp/n1s;
                look_up_n1[ai][si]=(ai-n1*si*1.0/M)*(ai-n1*si*1.0/M)/(n1*si*1.0/M);
                
            }            
        }        
    }
    cout <<"OK" <<endl <<endl;
    //cout <<n0 <<" " <<n1 <<" " <<M <<endl;    
    //for (int i=0 ; i<10 ; ++i)
    //    cout <<look_up_n0[0][i] <<endl;
    return 0;
}

//Use a look-up table to accelerate the computing efficiency. IN USE NOW.
double g_test(int a2, int a3, int A, int C, int S, int P, int G){    
	int a1,c1,c2,c3;	
	a1=A-a2-a3;
	c2=P-a2;
	c3=G-a3;	
	c1=C-c2-c3;
    
	double result=0;
    //cout <<"chi begin" <<endl;
    //cout <<"    " <<a1 <<" " <<a2 <<" " <<a3 <<" " <<c1 <<" " <<c2 <<" " <<c3 <<" " <<S <<" " <<P <<" " <<G <<" " <<endl;    
    result=look_up_n0[a1][S]+look_up_n0[a2][P]+look_up_n0[a3][G]+look_up_n1[c1][S]+look_up_n1[c2][P]+look_up_n1[c3][G];
    //cout <<"chi end" <<endl;
    
    //cout <<"    " <<a1 <<" " <<a2 <<" " <<a3 <<" " <<c1 <<" " <<c2 <<" " <<c3 <<" " <<S <<" " <<P <<" " <<G <" " <<result <<endl;    
    //cout <<result <<look_up_n0[a1][S] <<endl;    
    //cout <<"\t"<<result <<endl;
    //cout <<look_up_n0[a1][S] <<" " <<look_up_n0[a2][P] <<" " <<look_up_n0[a3][G] <<endl;
    //cout <<look_up_n1[c1][S] <<" " <<look_up_n1[c2][P] <<" " <<look_up_n1[c3][G] <<endl;
	return result;	
}


int create_g_test_lookup_tables(double **look_up_n0, double **look_up_n1)
{
    cout <<"Generating Lookup Tables...";
    
    int n1=NUM_OF_ONES_IN_Y;
    int n0=MAX_NUM_OF_INDIVIDUALS-n1;    
    int M=MAX_NUM_OF_INDIVIDUALS;
    double look_up_tmp;
    int n0s,n1s;    
    for (int si=0 ; si<MAX_NUM_OF_INDIVIDUALS+1 ; ++si){
        //n0s=n0*si;
        for (int ai=0 ; ai<MAX_NUM_OF_INDIVIDUALS+1 ; ++ai){
            if (ai==0){
                look_up_n0[ai][si]=0;                
            }
            else {                
	      //ook_up_tmp=M*ai-n0s;            
	      //look_up_n0[ai][si]=look_up_tmp*look_up_tmp/n0s;
	      look_up_n0[ai][si]=ai*log(ai*M*1.0/(si*n0));                
            }            
        }        
    }

    for (int si=0 ; si<MAX_NUM_OF_INDIVIDUALS+1 ; ++si){
        //n1s=n1*si;
        //look_up_n1[0][n1s]=0;        
        for (int ai=0 ; ai<MAX_NUM_OF_INDIVIDUALS+1 ; ++ai){
            if (ai==0)
                look_up_n1[ai][si]=0;
            else {                
	      //look_up_tmp=M*ai-n1s;            
	      //look_up_n1[ai][si]=look_up_tmp*look_up_tmp/n1s;
	      look_up_n1[ai][si]=ai*log(ai*M*1.0/(si*n1));
            }            
        }        
    }
    cout <<"OK" <<endl <<endl;
    //cout <<n0 <<" " <<n1 <<" " <<M <<endl;    
    //for (int i=0 ; i<10 ; ++i)
    //    cout <<look_up_n0[0][i] <<endl;
    return 0;
}
