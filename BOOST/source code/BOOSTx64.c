/* Copyright (C) (2010) (Can YANG) <eeyang@ust.hk>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

To collect contigency tables,
We first pre-calculate the number of 1's of 0~2^16 using the hamming weight to count the number of 1's in the string. (see http://en.wikipedia.org/wiki/Hamming_weight)
see the function bitCount. We store them in the vector "wordbits". Then use table looking method. see function popcount.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef long long   int64;
typedef unsigned long long uint64;
#define FMT_INT64   "%lld"
#define FMT_UINT64   "%llu"
#define FMT_HEX64   "%llx"

struct genotype{
	uint64 *genocase;
	uint64 *genoctrl;
};

struct MarginalDistr{
	int MarginalDistrSNP[3];
	int MarginalDistrSNP_Y[3][2];
};

static unsigned char wordbits[65536];// { bitcounts of ints between 0 and 65535 };
static int popcount( uint64 i )
{
    return( wordbits[i&0xFFFF] + wordbits[(i>>16)&0xFFFF] + wordbits[(i>>32)&0xFFFF] + wordbits[i>>48]);
}


double Abs(double a)
{
	return((a<0) ? -a : a);
}

int GetDataSize(char *filename, int **DataSize)
{
	FILE * fp, *fp_i;
	int c, ndataset;
	time_t st,ed;
	int n, p, i, flag,ii;
	char filename_i[100];


	fp = fopen(filename,"r");
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	ndataset = 0;
	while(!feof(fp)) {
		ndataset ++;
		fscanf(fp, "%s\n", &filename_i);
	}

	*DataSize = (int *)calloc( ndataset*2, sizeof(int));

	ii = 0;
	rewind(fp);
	while(!feof(fp)) {
		ii ++;
	fscanf(fp, "%s\n", &filename_i);

	fp_i = fopen(filename_i, "r");
	if(fp_i == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename_i);
		exit(1);
	}
	printf("start getting data size of file %d: %s\n", ii, filename_i);
	time(&st);
    //initialization
	if (ii == 1)
	{
	n = 0;//samples number

	// find the number of samples: n
	while(1)
	{
		int c = fgetc(fp_i);//read a character from the data file
		switch(c)
		{
		case '\n'://the end of line
			n++;
			break;
			// fall through,
			// count the '-1' element
		case EOF://file end
			goto out;
		default:
			;
		}
	}

	}
out:
	rewind(fp_i);//Repositions the file pointer to the beginning of a file

	// find number of variables: p
	p= 0;
	i= 0;
	flag = 1;
	while(1)
	{
		c = getc(fp_i);
		if(c=='\n') goto out2;//end of line
		if(isspace(c))
		{
			flag = 1;
		}
		/*do {
			c = getc(fp);
			if(c=='\n') goto out2;//end of line
		} while(isspace(c));//space
		*/
		if (!isspace(c) && (flag==1))
		{
			p++;//indicate the dimension of the vector
			flag = 0;
		}

	}
out2:
	fclose(fp_i);

	time(&ed);

//	DataSize[0] = n;
	(*DataSize)[ndataset * 0 + ii - 1] = n;
	(*DataSize)[ndataset * 1 + ii - 1] += p-1;

	}

	fclose(fp);
	//printf("Data contains %d rows and %d column. \n", n, p);

	printf("cputime for getting data size: %d seconds.\n", (int) ed - st);
	return ndataset;
}

//see http://en.wikipedia.org/wiki/Hamming_weight

int bitCount(uint64 i)
{
	i = i - ((i >> 1) & 0x5555555555555555);
	i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
	i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
	i = i + (i >> 8);
	i = i + (i >> 16);
	i = i + (i >> 32);
	return (int)i & 0x7f;
}

void CalculateMarginalEntropy(struct genotype *pgeno, int nsnp, int n, int nlongintcase, int nlongintctrl, double *MarginalEntropySNP, double *MarginalEntropySNP_Y)
{
    int i1, i2, i3;
    int count;
	double tmp, ptmp;
    int GenoMarginalDistr[3][2];

    for (i1 = 0; i1< nsnp; i1++)
    {
        for (i2 = 0; i2<3; i2++)
        {
            count = 0;
            for (i3 = 0; i3< nlongintcase; i3++)
			{
			    count += bitCount(pgeno[i1*3 + i2].genocase[i3]);
			}
			GenoMarginalDistr[i2][0] = count;

			count = 0;
			for (i3 = 0; i3< nlongintctrl; i3++)
			{
			    count += bitCount(pgeno[i1*3 + i2].genoctrl[i3]);
			}
			GenoMarginalDistr[i2][1] = count;
        }

        for (i2 = 0; i2<3; i2++)
        {
			tmp = (double) GenoMarginalDistr[i2][0] + GenoMarginalDistr[i2][1];
			if ( tmp > 0)
			{
				ptmp = tmp/n;
				MarginalEntropySNP[i1] += -(ptmp)*log(ptmp);
			}

			if (GenoMarginalDistr[i2][0]>0)
			{
				ptmp = (double) GenoMarginalDistr[i2][0]/n;
				MarginalEntropySNP_Y[i1] += -ptmp*log(ptmp);
			}

			if (GenoMarginalDistr[i2][1]>0)
			{
				ptmp = (double) GenoMarginalDistr[i2][1]/n;
				MarginalEntropySNP_Y[i1] += -ptmp*log(ptmp);
			}

        }
    }

}

void CalculateMarginalDistr(struct genotype *pgeno, int nsnp, int n, int nlongintcase, int nlongintctrl, struct MarginalDistr *pMarginalDistr)
{
    int i1, i2, i3;
    int count;
	//double tmp, ptmp;
    //int GenoMarginalDistr[3][2];

    for (i1 = 0; i1< nsnp; i1++)
    {
        for (i2 = 0; i2<3; i2++)
        {
            count = 0;
            for (i3 = 0; i3< nlongintcase; i3++)
			{
				count += bitCount(pgeno[i1*3 + i2].genocase[i3]);
			}
			pMarginalDistr[i1].MarginalDistrSNP_Y[i2][0] = count;

			count = 0;
			for (i3 = 0; i3< nlongintctrl; i3++)
			{
				count += bitCount(pgeno[i1*3 + i2].genoctrl[i3]);
			}
			pMarginalDistr[i1].MarginalDistrSNP_Y[i2][1] = count;

			pMarginalDistr[i1].MarginalDistrSNP[i2] = pMarginalDistr[i1].MarginalDistrSNP_Y[i2][0] + pMarginalDistr[i1].MarginalDistrSNP_Y[i2][1];
        }


    }

}


void CalculateGenoJointDistr(struct genotype *pgeno, int nsnp, int nLongIntcase, int nLongIntctrl, int *GenoDistr, int j1, int j2, struct MarginalDistr *pMarginalDistr)
{
	int i1, i2, i3;

	register int count;
//	uint64 tmp;

	for (i1 = 0; i1<2 ; i1++)
	{
		for (i2 = 0; i2 <2; i2++)
		{
		    count = 0;
			for (i3 = 0; i3< nLongIntcase; i3++)
			{
			    //count += bitCount(pgeno[j1*3 + i1].genocase[i3] & pgeno[j2*3 + i2].genocase[i3]);
				count += popcount(pgeno[j1*3 + i1].genocase[i3] & pgeno[j2*3 + i2].genocase[i3]);
				//GenoDistr[i1*3 + i2] += bitCount(pgeno[j1*3 + i1].genocase[i3] & pgeno[j2*3 + i2].genocase[i3]);
			}

            GenoDistr[i1*3 + i2] = count;
            count = 0;

			for (i3 = 0; i3< nLongIntctrl; i3++)
			{
			    //count += bitCount(pgeno[j1*3 + i1].genoctrl[i3] & pgeno[j2*3 + i2].genoctrl[i3]);
				count += popcount(pgeno[j1*3 + i1].genoctrl[i3] & pgeno[j2*3 + i2].genoctrl[i3]);
				//GenoDistr[9 + i1*3 + i2] += bitCount(pgeno[j1*3 + i1].genoctrl[i3] & pgeno[j2*3 + i2].genoctrl[i3]);
			}
			GenoDistr[9 + i1*3 + i2] = count;
		}

	}
//for case
	GenoDistr[2] = pMarginalDistr[j1].MarginalDistrSNP_Y[0][0] - GenoDistr[0] - GenoDistr[1];
	GenoDistr[5] = pMarginalDistr[j1].MarginalDistrSNP_Y[1][0] - GenoDistr[3] - GenoDistr[4];

	GenoDistr[6] = pMarginalDistr[j2].MarginalDistrSNP_Y[0][0] - GenoDistr[0] - GenoDistr[3];
	GenoDistr[7] = pMarginalDistr[j2].MarginalDistrSNP_Y[1][0] - GenoDistr[1] - GenoDistr[4];

	GenoDistr[8] = pMarginalDistr[j2].MarginalDistrSNP_Y[2][0] - GenoDistr[2] - GenoDistr[5];

//for ctrl
	GenoDistr[11] = pMarginalDistr[j1].MarginalDistrSNP_Y[0][1] - GenoDistr[9] - GenoDistr[10];
	GenoDistr[14] = pMarginalDistr[j1].MarginalDistrSNP_Y[1][1] - GenoDistr[12] - GenoDistr[13];

	GenoDistr[15] = pMarginalDistr[j2].MarginalDistrSNP_Y[0][1] - GenoDistr[9] - GenoDistr[12];
	GenoDistr[16] = pMarginalDistr[j2].MarginalDistrSNP_Y[1][1] - GenoDistr[10] - GenoDistr[13];

	GenoDistr[17] = pMarginalDistr[j2].MarginalDistrSNP_Y[2][1] - GenoDistr[11] - GenoDistr[14];

}



int main()
{
	time_t st, ed;

	struct genotype * pgeno;
	struct MarginalDistr *pMarginalDistr;

	/* Declare variable */
	int *DataSize;
	int ndataset;

    int n, p;  //n--number of samples; p number of varibles
	int ncase, nctrl, nlongintcase, nlongintctrl;
	int icase, ictrl;

	int *GenoJointDistr, *AlleleJointDistr;
	double or_aff, v_aff, or_unf, v_unf, pvalPLINK;//PLINK
	double *zval; //PLINK

    double *MarginalEntropySNP, *MarginalEntropySNP_Y, *MarginalAssociation;
	double JointEntropyTwoSNP, JointEntropyTwoSNP_Y, MarginalEntropyY, ptmp1, ptmp2;
	double InteractionMeasure;
	double *InteractionMeasureSNPpair;
	double maxInteraction, minInteraction;
	double *Pab, *Pbc, *Pca; // conditional probability P(a|b), P(b|c), P(c|a)
	double tao; // normalization term;


	int *DistrCollection;
	int *InteractionSNPpairs;
	int InteractionCount;

	int flag, i, ii, j, k, j1,j2;
	int c, tmp;
  // used for post-correction (post-correction is exact solution of loglinear model)
	static double mu[3][3][2] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	static double mu0[3][3][2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	static double mutmp[3][3][2];
	static double mu0tmp[3][3][2];

	static double mu_ij[3][3];
	static double mu_ik[3][2];
	static double mu_jk[3][2];
	double muError;
    
	int LengthLongType=64;
	uint64 mask1 = 0x0000000000000001;

	int buffersize = 50000;
	double thresholdRecord = 30;
	// load data
	FILE *fp, *fp_i;
	char filename_i[100];
	char filename[100] = "filenamelist.txt";
	fp = fopen(filename,"r");
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		return 1;
	}

	// precompute the wordbits (a global variable)
    for (i = 0; i<65536; i++)
    {
        wordbits[i] = bitCount(i);
        //printf("%d\n",wordbits[i]);
    }
	printf("start loading ...\n");

	ndataset = GetDataSize(filename, &DataSize);

	n = DataSize[0];
	p = 0;
	printf("n = %d\n", n);
	for (i = 0; i<ndataset ; i++)
	{
		p += DataSize[ndataset*1 + i];
		printf("DataSize %d-th file: p[%d] = %d \n", i+1, i+1, DataSize[ndataset*1 + i]);
	}
	printf("p = %d\n",p);


	// get ncase and nctrl
	i = 0;
	j = 0;

	ncase = 0;
	nctrl = 0;

	rewind(fp);

	// only use the first file to get ncase and nctrl
	fscanf(fp, "%s\n", &filename_i);
	printf("%s\n", filename_i);
	fp_i = fopen(filename_i, "r");

	while(!feof(fp_i)) {

		//if (n*j + i == 400)
		//{
		//	printf("%d,%d\n",i,j);
		//}

		/* loop through and store the numbers into the array */
		if(j==0)
		{
			//j = 0 means read ind class label y
			fscanf(fp_i, "%d", &tmp);

			if (tmp)
			{
				// tmp=1 means case
				ncase++;

			}
			else
			{
				nctrl ++;

			}
			j++;
		}
		else
		{
			fscanf(fp_i, "%d", &tmp);
			j++; //column index
			if (j==(DataSize[ndataset]+1)) // DataSize[ndataset] is the nsnp in the first dataset
			{
				j=0;
				i++; // row index
			}

		}

		if (i>=n)
		{
			break;
		}
	}

	printf("total sample: %d (ncase = %d; nctrl = %d).\n", n, (int)ncase, (int)nctrl);

	nlongintcase = ceil( ((double) ncase)/LengthLongType);
	nlongintctrl = ceil( ((double) nctrl)/LengthLongType);
	printf("nLongIntcase = %d; nLongIntctrl = %d.\n", nlongintcase, nlongintctrl);

	//calloc memory for bit representation
	pgeno = (struct genotype *)malloc(sizeof(struct genotype) * p * 3);// p SNPs, each contains 3 genotypes
	for (j = 0; j< 3*p ; j++)
	{
		(pgeno+j)->genocase = (uint64 *)calloc( nlongintcase, sizeof(uint64));
		(pgeno+j)->genoctrl = (uint64 *)calloc( nlongintctrl, sizeof(uint64));
	}

	//load data to bit representation

	
	rewind(fp);

	time(&st);
	j = 0; // column index
	ii = 0; // file index
	k = 0;
	while(!feof(fp)) { 
		ii ++;
		fscanf(fp, "%s\n", &filename_i);
		
		fp_i = fopen(filename_i, "r");
		if(fp_i == NULL)
		{
			fprintf(stderr,"can't open input file %s\n", filename_i);
			exit(1);
		}
		
		i = 0; //row index
		icase = -1;
		ictrl = -1;

		printf("Loading data in file %d: %s\n", ii, filename_i);
		while(!feof(fp_i)) { 
			/* loop through and store the numbers into the array */

			if(j==0)
			{
				//j = 0 means read class label y
				fscanf(fp_i, "%d", &tmp);
				
				if (tmp)
				{
					// tmp=1 means case
					icase ++;
					flag = 1;
				}
				else
				{
					ictrl ++;
					flag = 0;
				}
				j++;
			}
			else
			{
				fscanf(fp_i, "%d", &tmp);
				
				if (flag)
				{
					pgeno[(j+k-1)*3 + tmp].genocase[icase/LengthLongType] |= (mask1 << (icase%LengthLongType));
				}
				else
				{
					
					pgeno[(j+k-1)*3 + tmp].genoctrl[ictrl/LengthLongType] |= (mask1 << (ictrl%LengthLongType));
				}
				
				j++; //column index
				if (j==(DataSize[ndataset + ii-1]+1))
				{
					j=0;
					i++; // row index
				}
				
			}
			
			if (i>=n)
			{
				break;
			}
		}
		
		fclose(fp_i);
		k += DataSize[ndataset + ii-1];
	}
	
	fclose(fp);
	//printf("Number of numbers read: %d\n\n", n*p);
	time(&ed);
	printf("cputime for loading data: %d seconds\n", (int)ed -st);

	free(DataSize);


	// calculate marginal distribution

	ptmp1 = (double) ncase/n;
	MarginalEntropyY =  -ptmp1 *log(ptmp1) - (1-ptmp1) *log(1-ptmp1);

	MarginalEntropySNP = (double *)calloc(p, sizeof(double));
	MarginalEntropySNP_Y = (double *)calloc(p, sizeof(double));
	MarginalAssociation = (double *)calloc(p, sizeof(double));

    CalculateMarginalEntropy(pgeno, p, n, nlongintcase, nlongintctrl, MarginalEntropySNP, MarginalEntropySNP_Y);


    for (i = 0; i<p; i++)
    {
        MarginalAssociation[i] = (-MarginalEntropySNP_Y[i] + MarginalEntropySNP[i] + MarginalEntropyY)*n*2;
    }

    fp = fopen("MarginalAssoc.txt","w");
    for (i = 0; i < p; i++)
    {
		fprintf(fp,"%7d %f\n", i,MarginalAssociation[i]);
    }
    fclose(fp);

	pMarginalDistr = (struct MarginalDistr *)malloc(p *sizeof(struct MarginalDistr));
	CalculateMarginalDistr(pgeno, p, n, nlongintcase, nlongintctrl, pMarginalDistr);

    time(&st);
	// joint distribution
	GenoJointDistr = (int *)calloc(2*9, sizeof(int));
	AlleleJointDistr = (int *)calloc(2*4, sizeof(int) ); // for PLINK

	Pab		= (double *)calloc(9, sizeof(double));
	Pbc		= (double *)calloc(6, sizeof(double));
	Pca		= (double *)calloc(6, sizeof(double));

	DistrCollection = (int *)calloc(1001,sizeof(int));
	InteractionSNPpairs = (int *) calloc(buffersize*2, sizeof(int));
	InteractionMeasureSNPpair = (double *) calloc(buffersize,sizeof(double));

	InteractionCount = 0;

	maxInteraction = - 9999999;
	minInteraction = 9999999;


	for (j1 = 0; j1< p-1; j1++)
	//for (j1 = 0; j1< 5231; j1++)
	{
		for (j2 = j1+1; j2<p; j2++)
		//for (j2 = 5231; j2<9632; j2++)
		{

			// calculate joint distribution
			CalculateGenoJointDistr(pgeno, p, nlongintcase, nlongintctrl, GenoJointDistr, j1, j2, pMarginalDistr);

			/*
GenoJointDistr: the index is as follows:
					AABB
			Case	0	1	2	3	4	5	6	7	8
			Ctrl	9	10	11	12	13	14	15	16	17
*/
			// P(A|B)
			/* index
					B = 0   B = 1	B = 2
			A = 0	0		3		6
			A = 1	1		4		7
			A = 2	2		5		8
			*/

			// i for B, j for A
			for (i = 0; i< 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					Pab[3*i+j] = (double) (GenoJointDistr[3*j+i] + GenoJointDistr[3*j+i+9])/ pMarginalDistr[j2].MarginalDistrSNP[i];
				}

			}

			//P(B|C)
			/* index
					C = 0   C = 1
			B = 0	0		3
			B = 1	1		4
			B = 2	2		5
			*/
			//i for C (Y, class label), j for B

			i = 0;
				for (j = 0; j<3; j++)
				{
					Pbc[3*i + j] = (double) pMarginalDistr[j2].MarginalDistrSNP_Y[j][i]/ncase;
				}
			i = 1;
				for (j = 0; j<3; j++)
				{
					Pbc[3*i + j] = (double) pMarginalDistr[j2].MarginalDistrSNP_Y[j][i]/nctrl;
				}

			//P(C|A)
			/*index
					A=0		A=1		A=2
			C = 0    0		2		4
			C = 1	 1		3		5
			*/
			for (i = 0; i<3; i++)// i for A
			{
				for (j = 0; j<2; j++) // j for C
				{
					Pca[2*i + j] = (double) pMarginalDistr[j1].MarginalDistrSNP_Y[i][j]/ pMarginalDistr[j1].MarginalDistrSNP[i];
				}
			}

			//Papprx = Pab*Pbc*Pca
			// tao = sum(Papprx)
			// sum(p.* log(p) - p.*log(Pappr) + p.* log(tao))
			tao = 0.0;
			InteractionMeasure=0.0;
			for (i = 0; i<3; i++)// index for A
			{
				for (j = 0; j<3; j++) //index for B
				{
					for (k = 0; k<2; k++)	//index for C
					{
						ptmp1 = (double)GenoJointDistr[k*9 + i*3 +j]/n;
						if (ptmp1>0)
						{
							InteractionMeasure += ptmp1 * log(ptmp1);
						}

						ptmp2 = Pab[3*j+i]*Pbc[3*k+j]*Pca[2*i+k];
						if (ptmp2>0)
						{
							InteractionMeasure += -ptmp1*log(ptmp2);
							tao += ptmp2;
						}
					}
				}
			}

			InteractionMeasure = (InteractionMeasure+log(tao))*n*2;

			if (InteractionMeasure > maxInteraction)
			{
				maxInteraction = InteractionMeasure;
			}

			if (InteractionMeasure < minInteraction)
			{
				minInteraction = InteractionMeasure;
			}

			if (InteractionMeasure > thresholdRecord)
			{

				if (InteractionCount >= buffersize)// buffersize is not enough, calloc new memory
				{
					buffersize = buffersize * 2;
					InteractionSNPpairs = (int *)realloc(InteractionSNPpairs, buffersize*2*sizeof(int));
					InteractionMeasureSNPpair = (double *)realloc(InteractionMeasureSNPpair, buffersize*sizeof(double));

					InteractionSNPpairs[2*InteractionCount] = j1;
					InteractionSNPpairs[2*InteractionCount + 1] = j2;

					InteractionMeasureSNPpair[InteractionCount] = InteractionMeasure;
					InteractionCount ++;
					//printf("InteractionCount: %6d\t(%6d,%6d)\t %f\n", InteractionCount,j1,j2, InteractionMeasure);
				}
				else
				{
					InteractionSNPpairs[2*InteractionCount] = j1;
					InteractionSNPpairs[2*InteractionCount + 1] = j2;

					InteractionMeasureSNPpair[InteractionCount] = InteractionMeasure;
					InteractionCount ++;
					//printf("InteractionCount: %6d\t(%6d,%6d)\t %f\n", InteractionCount,j1,j2, InteractionMeasure);
				}
			}

			// collect distribution based on KSA
			if (InteractionMeasure > 100)
			{
				DistrCollection[1000] ++;
				//printf("pair:%d,%d\t %f\n", j1,j2, InteractionMeasure);
			}
			else if(InteractionMeasure >0)
			{
				DistrCollection[(int)(InteractionMeasure/0.1)] ++;
			}


		}
		if ((j1+1)%100==0)
            printf("iteration %d.\n", j1+1);
	}

	//post-correction

	printf("Start post-correction: Exact Gtest...\n");

	zval = (double *)calloc(InteractionCount, sizeof(double));

	for (ii = 0; ii < InteractionCount ; ii++)
	{
		j1 = InteractionSNPpairs[2*ii];
		j2 = InteractionSNPpairs[2*ii + 1];
		CalculateGenoJointDistr(pgeno, p-1, nlongintcase, nlongintctrl, GenoJointDistr, j1, j2, pMarginalDistr);

		memcpy(mutmp, mu, 18*sizeof(double));
		memcpy(mu0tmp, mu0, 18*sizeof(double));

		//Iterative Proportional Fitting homogeneous model (see section 8.7.2 of Categorical Data Analysis (2nd Ed.))
		muError = 0.0;
		for (i = 0; i<3; i++)
		{
			for(j = 0;j <3; j++)
			{
				for (k = 0; k <2; k++)
				{
					muError += Abs(mutmp[i][j][k]-mu0tmp[i][j][k]);
				}
			}
		}

		while (muError > 0.001)
		{
			memcpy(mu0tmp, mutmp, 18*sizeof(double)); //mu0tmp = mutmp;

			// mu_ij
			for (i = 0; i<3; i++)
			{
				for (j = 0; j<3; j++)
				{
					mu_ij[i][j] = mutmp[i][j][0] + mutmp[i][j][1];
				}

			}
			//mu_ijk = mu_ijk*n_ij/mu_ij
			for (i = 0; i<3; i++)
			{
				for(j = 0;j <3; j++)
				{
					for (k = 0; k <2; k++)
					{
						if (mu_ij[i][j]>0)
						{
							mutmp[i][j][k] = mutmp[i][j][k] * (GenoJointDistr[3*i + j]+GenoJointDistr[9 + 3*i + j])/mu_ij[i][j];
						}
						else
							mutmp[i][j][k] = 0;

					}
				}
			}

			// mu_ik
			for (i = 0; i<3; i++)
			{
				for (k = 0; k<2; k++)
				{
					mu_ik[i][k] = mutmp[i][0][k] + mutmp[i][1][k] + mutmp[i][2][k];
				}

			}
			//mu_ijk = mu_ijk*n_ik/mu_ik
			for (i = 0; i<3; i++)
			{
				for(j = 0;j <3; j++)
				{
					for (k = 0; k <2; k++)
					{
						// mu(i,j,k) = mu(i,j,k) * n_ik(i,k)/mu_ik(i,1,k);
						if (mu_ik[i][k] > 0)
						{
							mutmp[i][j][k] = mutmp[i][j][k] * (GenoJointDistr[9*k + 3*i ]+GenoJointDistr[9*k + 3*i + 1] +GenoJointDistr[9*k + 3*i + 2])/mu_ik[i][k];
						}
						else
							mutmp[i][j][k] = 0;
					}
				}
			}

			// mu_jk
			for (j = 0; j<3; j++)
			{
				for (k = 0; k<2; k++)
				{
					mu_jk[j][k] = mutmp[0][j][k] + mutmp[1][j][k] + mutmp[2][j][k];
				}

			}

			//mu_ijk = mu_ijk*n_jk/mu_jk
			for (i = 0; i<3; i++)
			{
				for(j = 0;j <3; j++)
				{
					for (k = 0; k <2; k++)
					{
						// mu(i,j,k) = mu(i,j,k) * n_jk(k,j)/mu_jk(1,j,k);
						if (mu_jk[j][k] > 0)
						{
							mutmp[i][j][k] = mutmp[i][j][k] * (GenoJointDistr[9*k + j]+GenoJointDistr[9*k + j + 3] +GenoJointDistr[9*k + j + 6])/mu_jk[j][k];
						}
						else
							mutmp[i][j][k] = 0;

					}
				}
			}
			//calculate Error
			muError = 0.0;
			for (i = 0; i<3; i++)
			{
				for(j = 0;j <3; j++)
				{
					for (k = 0; k <2; k++)
					{
						muError += Abs(mutmp[i][j][k]-mu0tmp[i][j][k]);
					}
				}
			}

		}// end for while

		tao = 0.0;
		InteractionMeasure=0.0;
		for (i = 0; i<3; i++)// index for A
		{
			for (j = 0; j<3; j++) //index for B
			{
				for (k = 0; k<2; k++)	//index for C
				{
					ptmp1 = (double)GenoJointDistr[k*9 + i*3 +j]/n;
					if (ptmp1>0)
					{
						InteractionMeasure += ptmp1 * log(ptmp1);
					}
					ptmp2 = mutmp[i][j][k]/n;
					if (ptmp2>0)
					{
						InteractionMeasure += -ptmp1*log(ptmp2);
						tao += ptmp2;
					}
				}
			}
		}

		InteractionMeasure = (InteractionMeasure+log(tao))*n*2;
		if (InteractionMeasure < 30)
		{
			//printf("SNPpair (%7d,%7d) : Approx: %f; Exact: %f\n", j1,j2, InteractionMeasureSNPpair[ii], InteractionMeasure);
		}

		InteractionMeasureSNPpair[ii] = InteractionMeasure; // update the interactionMeasure;


			//PLINK
	//     BB Bb  bb
//     AA  a  b  c
//     Aa  d  e  f
//     aa  g  h  i
//          B            b
//     A  4a+2b+2d+e   4c+2b+2f+e
//     a  4g+2h+2d+e   4i+2h+2f+e

//  case
    AlleleJointDistr[0] = 4*GenoJointDistr[0] + 2*GenoJointDistr[1] + 2*GenoJointDistr[3] + GenoJointDistr[4];
    AlleleJointDistr[1] = 4*GenoJointDistr[2] + 2*GenoJointDistr[1] + 2*GenoJointDistr[5] + GenoJointDistr[4];
    AlleleJointDistr[2] = 4*GenoJointDistr[6] + 2*GenoJointDistr[7] + 2*GenoJointDistr[3] + GenoJointDistr[4];
    AlleleJointDistr[3] = 4*GenoJointDistr[8] + 2*GenoJointDistr[7] + 2*GenoJointDistr[5] + GenoJointDistr[4];
// control
    AlleleJointDistr[4] = 4*GenoJointDistr[9] + 2*GenoJointDistr[10] + 2*GenoJointDistr[12] + GenoJointDistr[13];
    AlleleJointDistr[5] = 4*GenoJointDistr[11] + 2*GenoJointDistr[10] + 2*GenoJointDistr[14] + GenoJointDistr[13];
    AlleleJointDistr[6] = 4*GenoJointDistr[15] + 2*GenoJointDistr[16] + 2*GenoJointDistr[12] + GenoJointDistr[13];
    AlleleJointDistr[7] = 4*GenoJointDistr[17] + 2*GenoJointDistr[16] + 2*GenoJointDistr[14] + GenoJointDistr[13];
//
    or_aff = log( (double)(AlleleJointDistr[0]*AlleleJointDistr[3])/ (double)(AlleleJointDistr[1]*AlleleJointDistr[2]) );
    v_aff = 1/(double)AlleleJointDistr[0] + 1/(double)AlleleJointDistr[1] + 1/(double)AlleleJointDistr[2] + 1/(double)AlleleJointDistr[3];

    or_unf = log( (double)(AlleleJointDistr[4]*AlleleJointDistr[7])/ (double)(AlleleJointDistr[5]*AlleleJointDistr[6]) );
    v_unf = 1/(double)AlleleJointDistr[4] + 1/(double)AlleleJointDistr[5] + 1/(double)AlleleJointDistr[6] + 1/(double)AlleleJointDistr[7];

    //zval[ii] = Abs( (or_aff - or_unf) / sqrt ( v_aff + v_unf ) );
	zval[ii] = ( (or_aff - or_unf) / sqrt ( v_aff + v_unf ) );


		if ((ii+1)%10000==0)
			printf("iteration %d.\n", ii+1);

	}
   
	time(&ed);

	printf("maxInteraction : %f \t minInteraction: %f \n", maxInteraction, minInteraction);

/*	fp = fopen("DistrCollection.txt","w");
	for (i = 0; i<1001; i++)
	{
		fprintf(fp, "%d\n", DistrCollection[i]);
	}

	fclose(fp);
*/
    fp = fopen("InteractionRecords.txt","w");
    for (i = 0; i < InteractionCount; i++)
    {
    	if (InteractionMeasureSNPpair[i] > thresholdRecord)
        fprintf(fp,"%7d\t%7d\t%7d\t%f\t%f\t%f\t%f\n", i,InteractionSNPpairs[2*i],InteractionSNPpairs[2*i+1], MarginalAssociation[InteractionSNPpairs[2*i]], MarginalAssociation[InteractionSNPpairs[2*i+1]], InteractionMeasureSNPpair[i], zval[i]);
    }
    fclose(fp);

/*	CalculateGenoJointDistr(pgeno, p, nlongintcase, nlongintctrl, GenoJointDistr, 0, 1, pMarginalDistr);


	printf("cputime: %d\n", (int)ed - st);
	for (i = 0; i<2; i++)
	{
		for (j = 0; j<9; j++)
		{
			printf("%3d\t", GenoJointDistr[i*9 + j]);
		}
		printf("\n");
	}
*/
	for (j = 0; j< 3*p ; j++)
	{
		free(pgeno[j].genocase);
		free(pgeno[j].genoctrl);
	}

	free(Pab);
	free(Pbc);
	free(Pca);


	free(DistrCollection);
	free(InteractionSNPpairs);
	free(InteractionMeasureSNPpair);

	free(MarginalEntropySNP);
	free(MarginalEntropySNP_Y);
	free(MarginalAssociation);
	free(GenoJointDistr);
	free(pgeno);
	free(pMarginalDistr);

	free(AlleleJointDistr);//PLINK
	free(zval);//PLINK


	return 1;
}
