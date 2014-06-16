library(mbmdr)

args <- commandArgs(TRUE)

filename <- args[1]
filenamepath <- args[2]
dataset <- read.csv(filename, head = T, sep = ",")

fit <- mbmdr(y=dataset$Y,data=dataset[,3:301],order=2)

save(fit, file = filenamepath)

print(fit)

mbmdr.PermTest(fit, 100, sig.level=1)