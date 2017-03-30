standardize <- function(genotypefile, company='23andMe', snpref='hgdp.base.txt')
{
   if (company=='23andMe') {
      X<-read.table(genotypefile)
   }
   if (company=='ftdna') {
      X<-read.table(genotypefile, sep=',',skip=1)
   }
   if (company=='geno2') {
      COUNT<-0
      TOKEN <- c(as.matrix(read.table(genotypefile,skip=COUNT,nrows=1)[1]))
      print(TOKEN)
      while (TOKEN!='[Data]') {
         TOKEN <- c(as.matrix(read.table(genotypefile,skip=COUNT,nrows=1)[1]))
         COUNT <- COUNT+1
      }
      X<-read.table(genotypefile,sep=',',skip=COUNT)[,c(2,4,5)]
      rownames(X)<-X[,1]
      REF <- read.table(snpref)
      rownames(REF)<-REF[,1]
      COMMON<-intersect(rownames(X),rownames(REF))
      REF <- REF[COMMON,]
      X<-X[COMMON,]
      X<-cbind(REF,paste(X[,2],X[,3],sep=''))
   }
   if (company=='geno2new') {
      X<-read.table(genotypefile,sep=',',skip=1)[,c(1,3,4)]
      DUP <- which(duplicated(X[,1]))
      X<-X[setdiff(1:dim(X)[1],DUP),]
      rownames(X)<-X[,1]
      REF <- read.table(snpref)
      rownames(REF)<-REF[,1]
      COMMON<-intersect(rownames(X),rownames(REF))
      REF <- REF[COMMON,]
      X<-X[COMMON,]
      X<-cbind(REF,paste(X[,2],X[,3],sep=''))
   }

   write.table(X[order(X[,1]),],file='genotype.txt',quote=F,row.names=F,col.names=F)
}
