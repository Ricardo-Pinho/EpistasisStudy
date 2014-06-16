args <- commandArgs(TRUE)

chisq <- as.numeric(args[1])
df <- as.numeric(args[2])

val <- 1-pchisq(chisq,df)

print(val)

