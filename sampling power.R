setwd("/Users/jannawilloughby/GDrive/gray bats - alabama")
library(scales)

####power analysis for mark recapture####
popsize = 2000000
collects = recollects = c(seq(500, 5000, 500), seq(7500, 50000, 5000))

#output dataframe
OUT = NULL

#iterate over number of collected samples (marks and recaptures)
for(c in 1:length(collects)){
  collect = collects[c]
  recollect = recollects[c]
  Nst = NULL
  #repeat each marking effort 100 times
  for(i in 1:1000){
    #create initial population
    pop1 = data.frame(uid = seq(1, popsize, 1), capture = rep(0, popsize), recapture = rep(0, popsize))
    #initial marking
    pop1$capture[(sample(x=c(1:popsize), size=collect, replace=F))] = 1
    #recapture
    pop1$recapture[(sample(x=c(1:popsize), size=recollect, replace=F))] = 1
    #estimate pop1 size and add to list
    Nst = c(Nst, ((sum(pop1$capture)+1)*(collect+1))/((sum(pop1$recapture[pop1$capture==1])+1)))
  }
  #record mean/sd for pop1 size estimates
  writeout = c(collect, recollect, mean(Nst), sd(Nst), quantile(Nst, probs=0.975), quantile(Nst, probs=0.025))
  #add these values to others
  OUT = rbind(OUT, writeout)
}
colnames(OUT) = c("collect", "recollect", "NstM", "NstSD", "NstUL", "NstLL")
rownames(OUT) = seq(1,nrow(OUT), 1)
o = as.data.frame(OUT)
o$NstSE = 1.96*o$NstSD/sqrt(i)
o$totalgenos = o$collect + o$recollect

#plot data nicely
off = 100
plot(-10000, -10000, xlim=c(0, 100000), ylim=c(0,3000000), xlab="total number of samples genotyped per year", ylab="estimated population size")
segments(x0=0, x1=100000, y0=popsize, y1=popsize, lty=2, col="grey50")
lines(x=o$totalgenos, y=o$NstM, lty=1, col=alpha("firebrick3", 0.5), lwd=2)
points(x=o$totalgenos, y=o$NstM, pch=19, col=alpha("firebrick3", 0.5), cex=1.5)
segments(x0=c(o$totalgenos - off), x1=c(o$totalgenos + off), y0=o$NstLL, y1=o$NstLL, col=alpha("firebrick3", 0.8))
segments(x0=c(o$totalgenos - off), x1=c(o$totalgenos + off), y0=o$NstUL, y1=o$NstUL, col=alpha("firebrick3", 0.8))
segments(x0=o$totalgenos, x1=o$totalgenos, y0=o$NstLL, y1=o$NstUL, col=alpha("firebrick3", 0.8))


