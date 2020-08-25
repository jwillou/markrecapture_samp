setwd("~/Desktop/")
library(scales)

####power analysis for mark recapture####
popsize = 20000
collects = recollects = seq(50, 1200, 50)

#output dataframe
OUT = NULL

#iterate over number of collected samples (marks and recaptures)
for(c in 1:length(collects)){
  collect = collects[c]
  recollect = recollects[c]
  Nst = NULL
  #repeat each marking effort 100 times
  for(i in 1:100){
    #create initial herd
    herd = data.frame(uid = seq(1, popsize, 1), capture = rep(0, popsize), recapture = rep(0, popsize))
    #initial marking
    herd$capture[(sample(x=c(1:popsize), size=collect, replace=F))] = 1
    #recapture
    herd$recapture[(sample(x=c(1:popsize), size=collect, replace=F))] = 1
    #estimate herd size and add to list
    Nst = c(Nst, ((sum(herd$capture)+1)*(collect+1))/((sum(herd$recapture[herd$capture==1])+1)))
  }
  #record mean/sd for herd size estimates
  writeout = c(collect, recollect, mean(Nst), sd(Nst))
  #add these values to others
  OUT = rbind(OUT, writeout)
}
colnames(OUT) = c("collect", "recollect", "NstM", "NstSD")
rownames(OUT) = seq(1,nrow(OUT), 1)
o = as.data.frame(OUT)
o$NstSE = 1.96*o$NstSD/sqrt(i)

#plot data nicely
off = 10
plot(-100, -100, xlim=c(min(collects)*2, max(collects)*2), ylim=c(0,25000), xlab="total number of samples genotyped per year", ylab="estimated population size")
segments(x0=min(o$collect*2), x1=max(o$collect*2), y0=popsize, y1=popsize, lty=2, col="grey50")
polygon(x=c(500-50, 500+50, 500+50, 500-50), y=c(0,0,25000,25000), col=alpha("goldenrod2", 0.5), border=NA)
lines(x=o$collect*2, y=o$NstM, lty=1, col=alpha("firebrick3", 0.5), lwd=2)
points(x=o$collect*2, y=o$NstM, pch=19, col=alpha("firebrick3", 0.5), cex=1.5)
segments(x0=c(o$collect*2 - off), x1=c(o$collect*2 + off), y0=o$NstM-o$NstSE, y1=o$NstM-o$NstSE, col=alpha("firebrick3", 0.8))
segments(x0=c(o$collect*2 - off), x1=c(o$collect*2 + off), y0=o$NstM+o$NstSE, y1=o$NstM+o$NstSE, col=alpha("firebrick3", 0.8))
segments(x0=o$collect*2, x1=o$collect*2, y0=o$NstM-o$NstSE, y1=o$NstM+o$NstSE, col=alpha("firebrick3", 0.8))

####collars####
popsize = 15000
collars = seq(10, 150, 10)
mortality = c(0.05, 0.10, 0.20, 0.30)

#output dataframe
OUT = NULL

#iterate over mortaility proportions
for(m in 1:length(mortality)){
  #iterate over number of collars deployed
  for(c in 1:length(collars)){
    #repeat 100 times for each parameter combination
    t = NULL
    for(r in 1:100){
      pop = data.frame(uid = seq(1, popsize, 1), collar = rep(0, popsize), dead = rep(0, popsize))
      pop$collar[sample(c(1:nrow(pop)), collars[c], replace=F)] = 1
      pop$dead[sample(c(1:nrow(pop)), round(nrow(pop)*mortality[m], 0), replace=F)] = 1
      t = c(t, nrow(pop[pop$collar==1 & pop$dead==1,,drop=F])/collars[c])
    }
    writeout = c(mortality[m], collars[c], mean(t, na.rm=T), sd(t, na.rm=T))
    OUT = rbind(OUT, writeout)
  }
}
o = as.data.frame(OUT)
colnames(o) = c("mortality", "collars", "estmortM", "estmortSD")
rownames(o) = seq(1, nrow(o), 1)

#plot data nicely
off = 1
colors4 = c("firebrick3", "dodgerblue3", "darkgoldenrod3", "darkorchid3")
ylimu = c(0.2, 0.3, 0.4, 0.5)
for(m in 1:length(mortality)){
  pdf(file=paste("mort", mortality[m], ".pdf"), width=5, height=5)
  plot(-100, -100, xlim=c(10,150), ylim=c(0,ylimu[m]), xlab="number of collars", ylab="estimated mortality rate")
  #polygon(x=c(60-4, 60+4, 60+4, 60-4), y=c(0,0,ylimu[m],ylimu[m]), col=alpha("goldenrod2", 0.5), border=NA)
  t = o[o$mortality==mortality[m],]
  segments(x0=min(o$collars), x1=max(o$collars), y0=mortality[m], y1=mortality[m], lty=2, col="grey50")
  lines(x=t$collars,  y=t$estmortM, lty=1, col=alpha(colors4[m], 0.5), lwd=2)
  points(x=t$collars, y=t$estmortM, pch=19, col=alpha(colors4[m], 0.5), cex=1.5)
  segments(x0=c(t$collars - off), x1=c(t$collars + off), y0=t$estmortM-t$estmortSD, y1=t$estmortM-t$estmortSD, col=alpha(colors4[m], 0.8))
  segments(x0=c(t$collars - off), x1=c(t$collars + off), y0=t$estmortM+t$estmortSD, y1=t$estmortM+t$estmortSD, col=alpha(colors4[m], 0.8))
  segments(x0=t$collars, x1=t$collars, y0=t$estmortM-t$estmortSD, y1=t$estmortM+t$estmortSD, col=alpha(colors4[m], 0.8))
  dev.off()
}


