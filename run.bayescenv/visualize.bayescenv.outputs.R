#Devon DeRaad
#24 april 2020
#visualizing bayescenv outputs

#X chrom
#read in file with extra info
X.pvals<-read.csv("~/Desktop/anopheles.scripz/X.pvals.csv")

#read in trace file for hum
library(coda)
chain <- read.table("~/Downloads/X.hum.sel" , header = TRUE)
chain <- mcmc(chain , thin = 10)   #Adapt thin to its actual value (10 is the default)
heidel.diag(chain)             #To test for convergence
effectiveSize(chain)           #To compute effective sample size
autocorr.diag(chain)           #To look for auto-correlation
plot(chain)                    #To plot the "trace" and the posterior distribution
dev.off()

#read in trace file for temp
chain <- read.table("~/Downloads/X.temp.sel" , header = TRUE)
chain <- mcmc(chain , thin = 10)   #Adapt thin to its actual value (10 is the default)
heidel.diag(chain)             #To test for convergence
effectiveSize(chain)           #To compute effective sample size
autocorr.diag(chain)           #To look for auto-correlation
plot(chain)                    #To plot the "trace" and the posterior distribution
dev.off()

#read in trace file for precip
chain <- read.table("~/Downloads/X.precip.sel" , header = TRUE)
chain <- mcmc(chain , thin = 10)   #Adapt thin to its actual value (10 is the default)
heidel.diag(chain)             #To test for convergence
effectiveSize(chain)           #To compute effective sample size
autocorr.diag(chain)           #To look for auto-correlation
plot(chain)                    #To plot the "trace" and the posterior distribution
dev.off()



#read in output results files
x.baye.hum<-read.table("~/Desktop/anoph.3.march.2020/geste/X.hum_fst.txt")
x.baye.precip<-read.table("~/Desktop/anoph.3.march.2020/geste/X.precip_fst.txt")
x.baye.temp<-read.table("~/Desktop/anoph.3.march.2020/geste/X.temp_fst.txt")

#visualize
dev.off()
par(mfrow=c(2,1))
plot(X.pvals$pos, -log10(x.baye.hum$qval_g))
plot(X.pvals$pos, x.baye.hum$fst)
plot(X.pvals$pos, -log10(x.baye.temp$qval_g))
plot(X.pvals$pos, x.baye.temp$fst)
plot(X.pvals$pos, -log10(x.baye.precip$qval_g))
plot(X.pvals$pos, x.baye.precip$fst)




baye.hum<-read.table("~/Desktop/anoph.3.march.2020/")



