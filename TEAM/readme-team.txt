== TEAM == 

[How to compile]
  0) This program is running in Linux system. The following packages may be required:
     boost-1.35, cmake, g++-4
  1) Get into the TEAM directory, and run "cmake .".
  2) then run "make"

[Usage] 
  ./team.sh geno_fn pheno_fn #individuals #SNPs #permutations fde_threshold

[Options Details]
geno_fn: 	the genotype data file name (the rows are SNPs, the columns are individuals)
pheno_fn: 	the phenotype data file name
#individual: 	the number of individuals in the data
#SNPs: 		the number of SNPs in the data
#permutations:	the number of permutations used in the significant test
#fdr_threshold:	the FDR threshold for significance

[Example]
  ./team.sh ../data/genotypes_1k_2520.txt ../data/phenotypes_2520.txt 2520 1000 100 0.2

[Output Description]

  ======Significant SNP-pairs=====     <---- List of Significant SNP-pairs
  SNP-pair (284,70):42.5534            <---- SNP-pair (Xi,Xj): Test score
  SNP-pair (284,76):38.5984
  SNP-pair (284,74):41.0445
  SNP-pair (283,70):38.6654
  ================

Ver 0.0.4
@All Rights Reserved.
