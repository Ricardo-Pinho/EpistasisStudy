#! /bin/bash

#Change to the path of the script.
script_path="$(readlink -f $(dirname "$0"))"
#previous="$(pwd)"
#cd "$current"


if [ "$#" -lt "6" ]; then
   echo "Not enough parameters!";
   echo "Usage: $0 tri_geno bi_pheno #individual #SNPs #permutations #qvalue_threshold";
   exit 1;
fi

#Generate a random number as the intermedia file name
r=`dd if=/dev/urandom count=1 2> /dev/null | cksum | cut -f1 -d" "`
#echo $r;

${script_path}/test_all -if_geno $1 -if_pheno $2 -of_pvalue $r.txt.qvalue -n_inds $3 -n_snps $4 -n_perms $5
#./pvalue $r.txt.000 >$r.txt.qvalue
${script_path}/get_snps -if_geno $1 -if_pheno $2 -if_qvalue $r.txt.qvalue -n_inds $3 -n_snps $4 -qvalue $6
rm $r.txt.qvalue

echo "@@100";
echo "Done!";

