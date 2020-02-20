###R Script to calculate Spearman correlations
#Import data table
arid.soil<-read.table(file.choose(), sep="\t", header=T) 

#Create a for loop for testing multiple correlations
#TOC is being variable tested against different nitrogen forms 
for (i in 16:19) {
  corr.result<-cor.test(arid.soil[,i], arid.soil$TOC, method="spearman")
  cat(paste('\nDependent var:', colnames(arid.soil)[i], '\n'))
  print(corr.result)
}

##Correlation between gene abundance and different nitrogen forms
#Correlation between 16S gene abundance and different nitrogen forms
for (i in 16:19) {
  corr.result1<-cor.test(arid.soil[,i], arid.soil$X16S.logcopy, method="spearman")
  cat(paste('\nDependent var:', colnames(arid.soil)[i], '\n'))
  print(corr.result)
}

#Correlation between amoA gene abundance and different nitrogen forms
for (i in 16:19) {
  corr.result2<-cor.test(arid.soil[,i], arid.soil$amoA.logcopy, method="spearman")
  cat(paste('\nDependent var:', colnames(arid.soil)[i], '\n'))
  print(corr.result)
}

#Correlation between ureC gene abundance and different nitrogen forms
for (i in 16:19) {
  corr.result3<-cor.test(arid.soil[,i], arid.soil$ureC.logCopy, method="spearman")
  cat(paste('\nDependent var:', colnames(arid.soil)[i], '\n'))
  print(corr.result3)
}




