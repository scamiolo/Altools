### R code from vignette source 'vignettes/DNAcopy/inst/doc/DNAcopy.Rnw'

###################################################
### code chunk number 1: DNAcopy.Rnw:74-75
###################################################
library(DNAcopy)


###################################################
### code chunk number 2: DNAcopy.Rnw:78-79
###################################################
dataTable <- read.table("logRatioFile",header=TRUE)
dataFrame <- data.frame(dataTable)


###################################################
### code chunk number 3: DNAcopy.Rnw:85-88
###################################################
CNA.object <- CNA(dataFrame$logRatio,dataFrame$Chromosome,dataFrame$Position,data.type="logratio",sampleid = "testSample")


###################################################
### code chunk number 4: DNAcopy.Rnw:96-97
###################################################
smoothed.CNA.object <- smooth.CNA(CNA.object)


###################################################
### code chunk number 5: DNAcopy.Rnw:105-106
###################################################
segment.smoothed.CNA.object <- segment(smoothed.CNA.object,  verbose=1)


###################################################
### code chunk number 6: DNAcopy.Rnw:120-121
###################################################
plot(segment.smoothed.CNA.object, plot.type="w")


###################################################
### code chunk number 7: DNAcopy.Rnw:129-130
###################################################
plot(segment.smoothed.CNA.object, plot.type="c") 


###################################################
### code chunk number 8: DNAcopy.Rnw:157-158
###################################################
plot(segment.smoothed.CNA.object, plot.type="p")
write.table(segment.smoothed.CNA.object[2],file="cnv.txt")


###################################################
### code chunk number 9: DNAcopy.Rnw:169-172
###################################################
sdundo.CNA.object <- segment(smoothed.CNA.object, 
                             undo.splits="sdundo", 
                             undo.SD=3,verbose=1)
write.table(sdundo.CNA.object[2],file="cnv_smoothed.txt")


###################################################
### code chunk number 10: DNAcopy.Rnw:177-178
###################################################
plot(sdundo.CNA.object,plot.type="s")


