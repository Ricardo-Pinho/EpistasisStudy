int get_args(int argc, char * argv[]);

void usage();

void usage_get_snps();

int get_args_get_snps(int argc, char * argv[]);

double read_qvalue_list(ifstream &qvalue_infile, double qvalue_threshold);

double chi_square(int a2, int a3, int A, int C, int S, int P, int G);

int create_chi_square_lookup_tables(double **look_up_n0, double **look_up_n1);

double g_test(int a2, int a3, int A, int C, int S, int P, int G);

int create_g_test_lookup_tables(double **look_up_n0, double **look_up_n1);

//int read_genotypes(ifstream &geno_infile, vector<triSNP> &X);

//biPhenotype read_phenotypes(ifstream &pheno_infile);

void generate_permutation(int num_of_permutations, int num_of_ones, list <biPhenotype>& Ylist);

bool get_snps_partition(__int64 n_snps, __int64 n_partitions, int partition_id, 
                        __int64 & from_pos, __int64 & to_pos, __int64 & size);
