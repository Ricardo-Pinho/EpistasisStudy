source('screen_clean.R')

args <- commandArgs(TRUE)

genotype <- args[1]
phenotype <- args[2]
filenamepath <- args[3]

print(phenotype)

phenotypetable <- read.table(phenotype, header = TRUE)


pheno <- phenotypetable[,1]

geno <- read.table(genotype, header = TRUE)


sc1 = screenClean(geno, pheno, L = 300, K_pairs = 100, response = "binomial", standardize = TRUE, alpha =  0.05)

save(sc1, file = filenamepath)

print(sc1$snp_screen)

print(sc1$snp_screen2)

dim(sc1$clean)

sc1$clean

dim(sc1$final)

sc1$final

