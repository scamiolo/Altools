#plot coverage graph
table <- read.table("slidingFile",header=FALSE)
xLabel <- table[,1][1]
table <- read.table("slidingFile",header=TRUE)
plot(table[,1],table[,2],xlab=xLabel,ylab="Coverage in window")


#plot snp graph
plot(table[,1],table[,4],xlab=xLabel,ylab="SNP in window")

#plot indel graph
plot(table[,1],table[,6],xlab=xLabel,ylab="Indel in window")
