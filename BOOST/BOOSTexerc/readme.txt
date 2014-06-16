
This code runs on window XP 64bit.


INPUT

(1) filenamelist: specify the file name of the data. (For example, simdata.txt)

(2) the file name of the data. The first column is the class lable, taking its value in {0,1}. From the second column to
the last column, each column represents a SNP, taking vaules in {0,1,2}.

**please see our two examples for the input.

OUTPUT

(1)

InteractionRecords.txt

In our current released code, we only output those SNP pairs with ( InteractionBOOST>30 )

its format is as follows:

index	SNP1	SNP2	singlelocusAssoc1	singlelocusAssoc2	InteractionBOOST	InteractionPLINK


The index begins with 0.
The index of SNP begins with 0.


details:

Let x1 and x2 be the variable of SNP1 and SNP2, respectively. Let Y be the class label.

In R language, our code actually compute )

newdata = data.frame( x1 = factor(SNP1),x2 = factor(SNP2),y = classlabel )
fit0 = glm(y~x1, family= "binomial", data = newdata)
fit1 = glm(y~x2, family= "binomial", data = newdata)
fit01 = glm(y~x1+x2, family= "binomial", data = newdata)
fit2 = glm(y~x1+x2+x1*x2, family= "binomial", data = newdata)

singlelocusAssoc1 = fit0$null.dev-fit0$dev, df = 2
singlelocusAssoc2 = fit1$null.dev-fit1$dev, df = 2
InteractionBOOST  = fit01$dev-fit2$dev, df = 4

InteractionPLINK is the z value obtained by using the statistic of PLINK

In current version, we do not provide P-values but it can be easy to obtain using standard statistical software, e.g., R.
Then one can use either Bonferroni correction or FDR to control the false positive rate.


(2)
MarginalAssoc.txt

this file gives all single-locus association statistic value (df=2). This file has only two columns.

the format is 

SNPindex   single-locus test value

