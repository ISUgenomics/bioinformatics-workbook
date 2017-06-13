Select required columns for analysis. In this example, we are only interested in samples from B1 family. Make sure that te header line also has the same number of columns as the rest of the file.



```
awk '{print $1, $2, $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25}' all.counts.matrix > B1_family
awk '{print $1,$26, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49}' all.counts.matrix > B2_family
```


filein <-'B1_family'  ##Input file 
dataIn <- as.matrix(read.table(filein)) ##Read the data in file
dim(dataIn) ##Check if the number of columns is as expected
cols<-dim(dataIn)[2]
colsm1<-cols-1
dataIn2<-(dataIn[which(rowSums(dataIn[,1:cols])>colsm1),]) ##Retain only those rows that do not have low counts in all samples combined
dataIn3<-dataIn2[which(rowSums(sign(dataIn2[,1:cols]))>1),] ##Retain rows that have non-zero counts in at least one column 
dataIn<-as.matrix(dataIn3)
dim(dataIn) ##Number of transcripts that remained after clean up

##Quasiseq
upper.quartiles<-apply(dataIn,2,function(x) quantile(x,0.75)) ##Normalization - calculating upper quartile
scalingFactors<-mean(upper.quartiles)/upper.quartiles ##Generate scaling factor
range(scalingFactors)
dataIn.norm<-round((sweep(dataIn,2,scalingFactors,FUN="*"))) ##Generate normalized counts
trt <- as.factor(c(1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,2,2,2,2,2,2)) ##Assigning treatment groups to samples
mn <- as.numeric(rep(1,22)) ##Control
design.list<-vector('list',2)
design.list[[1]]<-model.matrix(~trt)
design.list[[2]]<-mn
log.offset<-log(apply(dataIn, 2, quantile,.75)) ##Normalization. calculate log offset
fit2<-QL.fit(round(dataIn), design.list,log.offset=log.offset, Model='NegBin') 
outfile-"B1_resvssus"
QLresultsPvaluePlot(fit2,outfile)

##EdgeR-QLFit
group <- factor(c(1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,2,2,2,2,2,2))
y <- DGEList(counts=dataIn,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tablrt <- topTags(qlf, n=Inf)
write.table(tablrt, file = "B1_resvssus_QLfit.txt")

##EdgeR-LRT
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
tablrt <- topTags(lrt, n=Inf)
write.table(tablrt, file = "B1_resvssus_LRT.txt")
</sxh>