WC_FST_FiniteSample_Haploids_2AllelesB_MCW<-function(AllCounts){
  #Input a matrix of the counts of each allele (columns) in each population (rows)
  #returns vector instead of list of Fst values, according to Weir
  
  n_pops<-dim(AllCounts)[1]
  r<-n_pops
  counts1 <- AllCounts[,1]
  sample_sizes <- rowSums(AllCounts)
  n_ave <- mean(as.numeric(sample_sizes))
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)  
  p_freqs = counts1/sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(p_freqs) in R
  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  #note: this differs slightly from var(p_freqs) in R
  
  T1 <- s2 - 1/(n_ave-1)*(p_ave*(1-p_ave) -(s2*(r-1)/r))
  T2 <- (n_c - 1)*(p_ave*(1-p_ave))/(n_ave-1) + (1 + (r-1)*(n_ave-n_c)/(n_ave-1))*s2/r
  
  FST <- T1/T2 
  
  return(c(He,p_ave, FST, T1, T2))
  
}

allelefreqchange.ind <- function(p0,p1, numind){
  p0.samp1 <- rbinom(1, numind, p0)/numind  
  p1.samp1 <- rbinom(1, numind, p1)/numind
  
  dp.true<- p0-p1
  dp.1<- p0.samp1-p1.samp1
  
  p.bar<- mean(p0,p1)
  p.bar.samp1<- mean(p0.samp1,p1.samp1)
  
  FST.0<-WC_FST_FiniteSample_Haploids_2AllelesB_MCW(
    matrix(c(p0*numind, p1*numind, (1-p0)*numind, (1-p1)*numind), ncol=2))
  
  FST.2<-WC_FST_FiniteSample_Haploids_2AllelesB_MCW(
    matrix(c(p0.samp1*numind, p1.samp1*numind, (1-p0.samp1)*numind, (1-p1.samp1)*numind), ncol=2))
    dp.true<- p0-p1
  dp.1<- p0.samp1-p1.samp1
  return(list(p0=p0,p1=p1,p0.samp1=p0.samp1,p1.samp1=p1.samp1, dp.true=dp.true, dp.1=dp.1, FST.0=FST.0[3], FST.2=FST.2[3], He.2=FST.2[1], T1=FST.2[4], T2= FST.2[5]))
}

allelefreqchange.pooled <- function(p0,p1, numind, numreads){
  p0.samp1 <- rbinom(1, numind, p0)/numind
  p0.samp2 <- rbinom(1, numreads, p0.samp1)/numreads
  
  p1.samp1 <- rbinom(1, numind, p1)/ numind
  p1.samp2 <- rbinom(1, numreads, p1.samp1)/numreads
  
  dp.true<- p0-p1
  dp.1<- p0.samp1-p1.samp1
  dp.2<- p0.samp2-p1.samp2 # this is observed allele frequency
  
  p.bar<- mean(p0,p1)
    p.bar.samp2<- mean(p0.samp2,p1.samp2)
  
  FST.0<-WC_FST_FiniteSample_Haploids_2AllelesB_MCW(
    matrix(c(p0, p1, (1-p0)*numreads, (1-p1)*numreads), ncol=2))
  
  FST.2<-WC_FST_FiniteSample_Haploids_2AllelesB_MCW(
    matrix(c(p0.samp2*numreads, p1.samp2*numreads, (1-p0.samp2)*numreads, (1-p1.samp2)*numreads), ncol=2))
  
  return(list(p0=p0,p1=p1, p0.samp2=p0.samp2, p1.samp2=p1.samp2, dp.true=dp.true, dp.1=dp.1, dp.2=dp.2, FST.0=FST.0[3], FST.2=FST.2[3], He.2=FST.2[1], T1=FST.2[4], T2= FST.2[5]))
}

l<-read.table("/Users/katie/Desktop/Rproj-OUTFLanK/MetaList_20130522_ED_v4_freqs.txt", header=TRUE)
head(l)
IM <- l[l$MODEL_SUBTYPE=="IM" & l$Neut10000==1,]
head(IM)
IM1 <- IM$X1_SSR1
IM2 <- IM$X2_SSR1
hist(IM$X1_SSR1)

### What I want to do here is compare the distribution of FST created by a IM simulations
### to the null distribution expected with the same sample size from sampling the same population
### without having undergone selection

  ### GET FST FROM SAME POP TWO TIME POINTS
      #quartz()
      p0<-sample(IM$X1_SSR1, 15000, replace=TRUE)
      p0<- runif(15000, 0.05,0.95)
      p1<-p0
      sample.size <- 46
      
      fp.23<- mapply(allelefreqchange.ind, p0,p1,sample.size)
      head(fp.23[,1:10],12)
      p0.samp1<- round(unlist(fp.23[3,])*sample.size,0)
      p1.samp1 <- round(unlist(fp.23[4,])*sample.size,0)
      FST.2times.overall<-sum(unlist(fp.23[10,]))/sum(unlist(fp.23[11,]))
      FST.2times<- unlist(fp.23[8,])
      He.2times<- unlist(fp.23[9,])
      #hist(FST.2times)
      #hist(p0.samp1)
      FST.2times.overall


  ### USE OVERALL FST TO IMPLEMENT IM SIMULATIONS
    setwd("~/Desktop/2014_AlleleFreqChangeProposal")
    #setwd("srcFDIST2inR")
    source("srcFDIST2inR/sourceAll.R")
    #LoadCFunctions_WithDelete(sourceWD) #uncomment this line to compile 1st time on your laptop
    loadCFunctions.withDelete("~/Desktop/2014_AlleleFreqChangeProposal/srcFDIST2inR")
    
    FST <-  FST.2times.overall
    number.reps <- 100000
    ndemes <- 2
    outfilepath.sims <- paste("IM_sims_Neut.txt", sep="")
    gens <- 1000 #get.num.gens(m,N)
        FST.set("CW93")
    IM.sims.replicate(ndemes, FST, sample.size, number.reps, outfilepath.sims, 
    diploid=FALSE, gens=gens, append.outfile=FALSE, progress=TRUE)
  FST.sim <- read.table(outfilepath.sims, header=TRUE)
  head(FST.sim)
  sum(FST.sim$numer)/sum(FST.sim$He)
  hist(FST.sim$FST)  
# write example dataset
# df<- data.frame(Ind=1:42*2, Pop=c(rep("Before", 42), rep("After", 42)), matrix(NA,ncol=15000,nrow=42*2))
# for (i in 1:15000){
#   s1<-c(rep(0, p0.samp1[i]), rep(1,(42-p0.samp1[i])))
#   s2 <- sample(s1,42,replace=FALSE)
#   df[1:42,(i+2)]<-s2
#   s3<- c(rep(0, p1.samp1[i]), rep(1,(42-p1.samp1[i])))
#   s4 <- sample(s3,42,replace=FALSE)
#   df[43:(42*2),(i+2)]<-s2  
# }
# write.table(df, "ExampleDataset42haploidsBeforeAfterSelection.txt", col.names=TRUE, row.names=FALSE)

### Now compare the above scenario (each individual sequenced) to the same 
### scenario with pooling
fp.46.40<-mapply(allelefreqchange.pooled, p0,p1,numind=sample.size,numreads=40)
head(fp.46.40[,1:10], 12)


### Compare plot: Almost perfectly the same for individuals, outliers for !
    br<- seq(-0.06, 0.5,0.01)
    hist.FST.2times<-hist(FST.2times, breaks=br, plot=FALSE)
    hist.FST.2pops<-hist(FST.sim$FST, breaks=br, plot=FALSE)
    hist.FST.2times.pooled<-hist(as.numeric(fp.46.40[9,]), breaks=br, plot=FALSE)

par(mfrow=c(2,1), mar=c(2,4,1,1), oma=c(2,2,0,0))   
    plot(hist.FST.2times$mids, sqrt(hist.FST.2times$dens), type="h", ylim=c(0,7), 
         yaxs="i", lwd=4, lend=2, bty="n", ylab="Sqrt(Density)", xlab="", col=rgb(0,0,1,0.5))
    points(hist.FST.2pops$mids+0.005, sqrt(hist.FST.2pops$dens), type="h", ylim=c(0,7), 
         yaxs="i", lwd=4, lend=2, ylab="Sqrt(Density)", bty="n", xlab=expression(F[ST]))
    legend(0.1,6.5, c("2 deme IM simulation", "Sequence 46 individuals separately" 
                     ), bty="n", 
           fill=c("black", rgb(0,0,1,0.5))) 

    plot(hist.FST.2times.pooled$mids, sqrt(hist.FST.2times.pooled$dens), type="h", ylim=c(0,7), 
         yaxs="i", lwd=4, lend=2, bty="n", ylab="Sqrt(Density)", xlab="", col=rgb(0,1,0,0.9))
    points(hist.FST.2pops$mids+0.005, sqrt(hist.FST.2pops$dens), type="h", ylim=c(0,7), 
         yaxs="i", lwd=4, lend=2, ylab="Sqrt(Density)", bty="n", xlab=expression(F[ST])) 
    legend(0.1,6.5, c("2 deme IM simulation", "Sequence 46 pooled individuals"), bty="n", 
           fill=c("black",  rgb(0,1,0,0.9))) 

   mtext(expression(F[ST]), outer=TRUE, side=1, line=1)

#library(hexbin)
#plot(hexbin(unlist(fp.23[9,]),unlist(fp.23[8,])))

h = seq(0,0.5,0.04)
CI95<-h
for (i in 2:length(h)){
  s1 <- FST[He < h[i] &  He > h[i-1]]
  s2 <- round(0.95*length(s1),0)
  CI95[i]<- sort(s1)[s2]
}
plot(jitter(unlist(fp.23[9,]), amount=0.02),jitter(unlist(fp.23[8,]),amount=0.02), pch=19,cex=0.5)
points(h,CI95, type="l", col="red", lwd=3)



fp.40.40<-mapply(allelefreqchange, p0,p1,40,40)
head(fp.40.40[,1:10], 10)
fp.40.1000<-mapply(allelefreqchange, p0,p1,40,1000)
fp.1000.40<-mapply(allelefreqchange, p0,p1,1000,40)

# fp.40.40<- replicate(1e04, allelefreqchange(0.5,0.5, 40, 40))
# dim(fp.40.40)
# head(fp.40.40[,1:10], 10)
# 
# fp.1000.40<- replicate(1e04, allelefreqchange(0.5,0.5, 1000, 40))
# fp.40.1000<- replicate(1e04, allelefreqchange(0.5,0.5, 40, 1000))

makeplot <- function(a){
  FST.samp<- unlist(a[9,])
  He.samp<- unlist(a[10,])
  hist(unlist(a[7,]), xlim=c(-0.5,0.5), breaks=seq(-0.6,0.6,0.05))
  #y<-hist(FST.samp, breaks=seq(-0.1,0.5,0.01), plot=FALSE)
  #plot(y$mids,sqrt(y$dens), type="h", xlim=c(-0.05,0.6), lwd=10, bty="n")
  plot(He.samp, FST.samp, xlim=c(0,0.5), ylim=c(0,0.7))
  head(FST.samp)
  #which(FST.samp==max(FST.samp))
}



quartz()
par(mfrow=c(3,2), mar=c(2,2,1,1))
makeplot(fp.40.40)
makeplot(fp.40.1000)
makeplot(fp.1000.40)

### NOTE:  IM SIMULATIONS BASED ON NUMBER OF INDIVIDUALS SAMPLED
### WOULD GIVE WRONG NULL DISTRIBUTION

par(mfrow=c(3,1))
hist(fp.40.40, xlim=c(-0.7,0.7), main="40 individuals and 40X coverage")
hist(fp.40.1000, xlim=c(-0.7,0.7), main="40 individuals and 1000X coverage")
hist(fp.1000.40, xlim=c(-0.7,0.7), main="1000 individuals and 40X coverage")

### Sig based on permutations

sample(fp.40.40[3,])




### Sig based on G-test
#G<- 2*sum(obs*ln(obs/exp))
reads=40
obs.p<-unlist(fp.40.40[4,])*reads
obs.q<-reads-obs.p
exp.p <- unlist(fp.40.40[3,])*reads
exp.q <- reads-exp.p
   
G<- 2*(obs.p*log(obs.p/exp.p) + obs.q*log(obs.q/exp.q))
hist(G)
likelihood.ratio.test

install.packages("Deducer", dependencies=TRUE)
library(Deducer)

likelihood.test(matrix(c(10,30,20,20), ncol=2))
?likelihood.test
