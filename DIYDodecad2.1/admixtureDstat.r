admixtureDstat <- function(parfile, listfile=c(),outgroup=c(),k=500, printzmatrix=F, outputfile='dstat.txt')
{
   PFILE <- c(as.matrix(read.table(parfile,skip=5,nrows=1)[1])) 
   AFILE <- c(as.matrix(read.table(parfile,skip=6,nrows=1)[1]))
   GFILE <- c(as.matrix(read.table(parfile,skip=2,nrows=1)[1]))
   NAMEFILE <- c(as.matrix(read.table(parfile,skip=4,nrows=1)[1]))
   PMATRIX <- read.table(PFILE)
   AMATRIX <- read.table(AFILE,row.names=1)
   GMATRIX <- read.table(GFILE,row.names=1)
   ELEMENTS <- c(as.matrix(read.table(NAMEFILE)))
   NELEMENTS <- length(ELEMENTS)
   COMMON <- intersect(rownames(AMATRIX),rownames(GMATRIX))
   print(paste(length(COMMON)," SNPs common between your genotype file ",GFILE," and the calculator file ",parfile))
   rownames(PMATRIX)<-rownames(AMATRIX)
   PMATRIX <- PMATRIX[COMMON,]
   AMATRIX <- AMATRIX[COMMON,]
   GMATRIX <- GMATRIX[COMMON,]
   GOOD <- which(c(as.matrix(GMATRIX[,3]))!='--')
   MISSING <- length(COMMON)-length(GOOD)
   print(paste(MISSING,' missing values in ',GFILE))   
   PMATRIX <- PMATRIX[GOOD,]
   AMATRIX <- AMATRIX[GOOD,]
   GMATRIX <- GMATRIX[GOOD,]
   print(paste("Analysis will be performed on ",dim(GMATRIX)[1]," SNPs"))
   A1 <- substr(c(as.matrix(GMATRIX[,3])),1,1)
   A2 <- substr(c(as.matrix(GMATRIX[,3])),2,2)
   MINOR <- c(as.matrix(AMATRIX[,1]))
   MAJOR <- c(as.matrix(AMATRIX[,2]))
   FLIPPED <- which((A1!=MINOR & A1!=MAJOR) & A1==A2)
   NOT_FLIPPED <- setdiff(1:dim(GMATRIX)[1],FLIPPED)
   GFREQ <- vector(length=dim(GMATRIX)[1])
   for (i in FLIPPED) {
      if (A1[i]=='A') {
         A1[i]<-'T'
         A2[i]<-'T'
      }
      else if (A1[i]=='T') {
         A1[i]<-'A'
         A2[i]<-'A'
      }
      else if (A1[i]=='C') {
         A1[i]<-'G'
         A2[i]<-'G'
      }
      else if (A1[i]=='G') {
         A1[i]<-'C'
         A2[i]<-'C'
      }
   }
   GFREQ[A1==MINOR & A1==A2] <- 0
   GFREQ[A1==MAJOR & A1==A2] <- 1
   GFREQ[A1!=A2]<-0.5

   print(paste(length(FLIPPED)," SNPs are flipped"))
   print(paste(sum(GFREQ==0.5)," SNPs are heterozygous"))


   ORDER <- order(GMATRIX[,1],GMATRIX[,2])
   PMATRIX <- PMATRIX[ORDER,]
   GFREQ <- GFREQ[ORDER]

   if (length(listfile)>0) {
      DSTATLIST<-read.table(listfile)
   }
   else {
      DSTATLIST <- array(dim=c((NELEMENTS-1)*(NELEMENTS-2),3))
      DSTATLIST[,3]<-outgroup
      WHICH_OUTGROUP <- which(ELEMENTS==outgroup)
      ITER <- 0
      for (i in setdiff(1:NELEMENTS,WHICH_OUTGROUP)) {
         for (j in setdiff(1:NELEMENTS,c(i,WHICH_OUTGROUP))) {
            ITER <- ITER+1
            DSTATLIST[ITER,1]<-ELEMENTS[i]
            DSTATLIST[ITER,2]<-ELEMENTS[j]
         }
      }
   }
   ZMATRIX <- array(dim=c(NELEMENTS,NELEMENTS))
   rownames(ZMATRIX)<-ELEMENTS
   colnames(ZMATRIX)<-ELEMENTS
   ZMATRIX[1:NELEMENTS,1:NELEMENTS] <- NA
   NDSTAT<-dim(DSTATLIST)[1]
   RESULT <- rbind(c('Pop1','Pop3','Outgroup','Dstat','Z'))
   for (COUNT in 1:NDSTAT) {
      if (NDSTAT>1) {
         pop1 <- which(c(as.matrix(DSTATLIST[COUNT,1]))==ELEMENTS)
         pop3 <- which(c(as.matrix(DSTATLIST[COUNT,2]))==ELEMENTS)
         pop4 <- which(c(as.matrix(DSTATLIST[COUNT,3]))==ELEMENTS)
      }
      else {
         pop1 <- which(c(as.matrix(DSTATLIST[1]))==ELEMENTS)
         pop3 <- which(c(as.matrix(DSTATLIST[2]))==ELEMENTS)
         pop4 <- which(c(as.matrix(DSTATLIST[3]))==ELEMENTS)
      }
      NUM <- (PMATRIX[,pop1]-GFREQ)*(PMATRIX[,pop3]-PMATRIX[,pop4])
      DENOM <- (PMATRIX[,pop1]+GFREQ-2*GFREQ*PMATRIX[,pop1])*(PMATRIX[,pop3]+PMATRIX[,pop4]-2*PMATRIX[,pop3]*PMATRIX[,pop4])
      DSTAT_GLOBAL <- sum(NUM)/sum(DENOM)
      NBLOCKS <- ceiling(length(GFREQ)/k)
      DSTAT <- vector(length=NBLOCKS)
      BLCK_LEN <- vector(length=NBLOCKS)
      BLCK_START <- 1
      ALL_SNPS <- 1:length(GFREQ)
      for (i in 1:NBLOCKS) {
          BLCK_END <- min(BLCK_START+k, length(GFREQ)) 
          BLCK_LEN[i] <- length(BLCK_START:BLCK_END)
          BLCK_SNPS <- setdiff(ALL_SNPS,BLCK_START:BLCK_END)
          NUM <- (PMATRIX[BLCK_SNPS,pop1]-GFREQ[BLCK_SNPS])*(PMATRIX[BLCK_SNPS,pop3]-PMATRIX[BLCK_SNPS,pop4])
          DENOM <- (PMATRIX[BLCK_SNPS,pop1]+GFREQ[BLCK_SNPS]-2*GFREQ[BLCK_SNPS]*PMATRIX[BLCK_SNPS,pop1])*(PMATRIX[BLCK_SNPS,pop3]+PMATRIX[BLCK_SNPS,pop4]-2*PMATRIX[BLCK_SNPS,pop3]*PMATRIX[BLCK_SNPS,pop4])
          DSTAT[i] <- sum(NUM)/sum(DENOM)
          BLCK_START<-BLCK_END+1
      }
      DSTAT_JACK <- NBLOCKS*DSTAT_GLOBAL-sum(DSTAT*(1-BLCK_LEN/length(GFREQ)))
      HJ <- length(GFREQ)/BLCK_LEN
      TAU <- HJ*DSTAT_GLOBAL-(HJ-1)*DSTAT
      VAR <- mean((TAU-DSTAT_JACK)^2/(HJ-1))
      RESULT <- rbind(RESULT,c(ELEMENTS[pop1],ELEMENTS[pop3],ELEMENTS[pop4], round(DSTAT_JACK,5), round(DSTAT_JACK/sqrt(VAR),2)))
      print(RESULT)
      if (printzmatrix) {
         ZMATRIX[pop1,pop3]<- round(DSTAT_JACK/sqrt(VAR),2)
         print(ZMATRIX)
      }
   }
   write.table(RESULT,file=outputfile,quote=F,row.names=F,col.names=F)
}
