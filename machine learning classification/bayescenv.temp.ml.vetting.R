#May 17 2020
#plot results from phase1 and verify with phase2 allele freq patterns
#devtools::install_github('exaexa/scattermore')
library(scattermore)
library(gridExtra)
library(VennDiagram)
library(qvalue)
library(outliers)
library(e1071)

all.sigs<-read.csv(file="~/Desktop/anoph.phase2/all.sigs.fdr.adj.csv")

#bring in phase1 allele frequencies
phase1.all.freqs<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.all.freqs.csv")
#make column for id
phase1.all.freqs$id<-paste(phase1.all.freqs$chrom,phase1.all.freqs$POS)

#bring in phase2 allele frequencies
phase2.all.freqs<-read.csv(file = "~/Desktop/anoph.phase2/phase2.all.freqs.csv")
#make id column
phase2.all.freqs$id<-paste(phase2.all.freqs$chrom,phase2.all.freqs$POS)

#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c(16,17,18,20,22)]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase1 <-metapops[order(metapops$latitude),]

#bring in dataframe with locality and environmental data for each individual phase2
phase2.alldata<-read.csv(file = "~/Downloads/phase2_localities_all_vairables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase2.alldata[,c(2,3,4,6,8)]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase2 <-metapops[order(metapops$latitude),]

#pull all significant bayescenv temp outliers
bayescenv.temp.outliers<- all.sigs[all.sigs$temp.q < .01,]
nrow(bayescenv.temp.outliers) #check number of outliers

#generate matching id's for the outlier dataset in question
sig.pos<-paste(bayescenv.temp.outliers$chrom,bayescenv.temp.outliers$pos)
#subset all AF values to only the outlier dataset, to save time rather than searching the entire 3M snp dataset
phase1.all.freqs.bayescenv.temp <- phase1.all.freqs[phase1.all.freqs$id %in% sig.pos,]

#build correlation dataframe for phase1 data
cor.frame<-data.frame(matrix(, nrow=0, ncol=10)) #init df
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase1.all.freqs.bayescenv.temp[phase1.all.freqs.bayescenv.temp$id == i,3:16]))
  mod<-summary(lm(frq ~ sort.pops.phase1$tmean)) #perform linear reg
  # store 9 values describing the correlation to feed into machine learning
  cor.frame<-rbind(cor.frame, cbind(i,
                                    max(abs(mod$residuals)),
                                    median(mod$residuals),
                                    mod$coefficients[1,1],
                                    mod$sigma,
                                    mod$r.squared,
                                    mod$fstatistic[1],
                                    cov(frq, sort.pops.phase1$tmean),
                                    cor(frq, sort.pops.phase1$tmean),
                                    sort(abs(frq-mean(frq)))[14]-sort(abs(frq-mean(frq)))[13]
  )
  )
  setTxtProgressBar(pb, j)
  j<-j+1
}

#fix row/colnames
colnames(cor.frame)<-c("id","maxresid","medresid","intercept","stdevresid","rsquared","f","cov","cor","outlier")
rownames(cor.frame) <- c()
#make numeric columns numeric
for (i in 2:10){
  cor.frame[,i]<-as.numeric(as.character(cor.frame[,i]))
  
}
#check
class(cor.frame$id)
class(cor.frame$medresid)

#calc correlation between these variables
round(cor(cor.frame[,-1]), 2)

#going to drop f, cov, and cor because they are highly correlated > .9
cor.frame<-cor.frame[,c(1:6,10)]

#pull out our random subset to hand identify as training dataset
cor.training<-cor.frame[sample(nrow(cor.frame),100),]
#read.csv(file="~/Downloads/bayescenv.temp.cor.training.csv")

par(mfrow=c(3,3))
#plot
sig.pos<-as.vector(cor.training$id)
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == i,3:16]))
  plot(sort.pops.phase1$tmean,frq, main = i, ylim =c(0,1))
  abline(lm(frq~sort.pops.phase1$tmean))
  setTxtProgressBar(pb, j)
  j<-j+1
}

#hand classify this training set
#"environment"=0, "outlier"=1
vec<-c(0,0,0,1,1,1,0,0,0,
       0,0,0,0,1,0,0,0,0,
       0,0,1,0,0,0,0,0,0,
       0,1,1,0,0,0,0,0,1,
       0,1,0,0,0,0,1,1,0,
       1,1,1,0,0,0,0,0,0,
       0,1,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,
       1,0,0,0,0,0,1,0,0,
       0,1,0,1,0,0,0,0,0,
       0,0,0,0,0,1,0,0,0,
       0)

vec[vec == 0]<-"environment"
vec[vec == 1]<-"outlier"

#add hand classification to the training set
cor.training$class<-as.factor(vec)
#write.csv(cor.training, file="~/Downloads/bayescenv.temp.cor.training.csv")

#train the model (not including max residual which decreased prediction accuracy)
m <- naiveBayes(cor.training[,2:7], cor.training[,8])
#check predictions against training data
table(predict(m, cor.training[,2:7]), cor.training[,8])
cor.training$pred.class<-predict(m, cor.training[,2:7]) #add prediction to df
cor.training$prob.env<-predict(m, cor.training[,2:7], type="raw")[,1] #add prediction probability to df
cor.training$prob.out<-predict(m, cor.training[,2:7], type="raw")[,2] #add prediction probability to df

#do prediction
cor.frame$class<-predict(m, cor.frame[,2:7]) #add prediction to df
cor.frame$prob.env<-predict(m, cor.frame[,2:7], type="raw")[,1] #add prediction probability to df
cor.frame$prob.out<-predict(m, cor.frame[,2:7], type="raw")[,2] #add prediction probability to df

#visualize prediction
par(mfrow=c(1,2))
pc.df<-as.data.frame(prcomp(cor.training[,2:7])$x)
plot(pc.df$PC1,pc.df$PC2, col=cor.training$class, main = "hand classified")
plot(pc.df$PC1,pc.df$PC2, col=predict(m, cor.training[,2:7]), main = "ML predicted")

#histograms of individual variables by predicted class
hist(cor.frame$intercept[cor.frame$class=="environment"], main="int",
     col = rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue"), ylim=c(0,1750))
hist(cor.frame$intercept[cor.frame$class=="outlier"],
     col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), add =TRUE)
hist(cor.frame$outlier[cor.frame$class=="environment"], main="outlier",
     col = rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue"), xlim=c(0,1))
hist(cor.frame$outlier[cor.frame$class=="outlier"],
     col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), add =TRUE)

#plot
plot(x=cor.frame$prob.env, y=cor.frame$intercept)
plot(x=cor.frame$prob.env, y=cor.frame$outlier)

#pca of all predicted points
par(mfrow=c(1,1))
pc.df<-as.data.frame(prcomp(cor.frame[,2:7])$x)
plot(pc.df$PC1,pc.df$PC2, col=cor.frame$class, main = "ML predicted")

#visualize specific training points
par(mfrow=c(3,3))
#plot
sig.pos<-as.vector(cor.training$id)
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == i,3:16]))
  if (cor.training$class[j] == cor.training$pred.class[j]){
    plot(sort.pops.phase1$tmean,frq, main = i, ylim =c(0,1),
         xlab=paste("prob env",round(cor.training$prob.env[j],2)), ylab=paste("prob out",round(cor.training$prob.out[j],2)))
    abline(lm(frq~sort.pops.phase1$tmean))
  }
  else{
    plot(sort.pops.phase1$tmean,frq, main = i, ylim =c(0,1), col = "red",
         xlab=paste("prob env",round(cor.training$prob.env[j],2)), ylab=paste("prob out",round(cor.training$prob.out[j],2)))
    abline(lm(frq~sort.pops.phase1$tmean), col="red")
  }
  setTxtProgressBar(pb, j)
  j<-j+1
}


#check prediction
table(cor.frame$class)
par(mfrow=c(1,1))
hist(cor.frame$prob.env)
#split confidently identified SNPs by class
cor.frame.outliers<-cor.frame[cor.frame$prob.out >.8,]
cor.frame.env<-cor.frame[cor.frame$prob.env >.8,]


#####step2 validate each class with phase2 data
####
###
##
#
#start with environmentally correlated SNPs
#generate matching id's for the outlier dataset in question
sig.pos<-as.vector(cor.frame.env$id)
#subset all AF values to only the outlier dataset, to save time rather than searching the entire 3M snp dataset
phase2.all.freqs.bayescenv.temp <- phase2.all.freqs[phase2.all.freqs$id %in% sig.pos,]
#generate list of SNP IDs found in phase2
sig.pos<-phase2.all.freqs.bayescenv.temp$id

#calc phase2 cor.frame
cor.frame.env.2<-data.frame(matrix(, nrow=0, ncol=10)) #init df
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase2.all.freqs.bayescenv.temp[phase2.all.freqs.bayescenv.temp$id == i,3:23]))
  mod<-summary(lm(frq ~ sort.pops.phase2$tmean)) #perform linear reg
  # store 9 values describing the correlation to feed into machine learning
  cor.frame.env.2<-rbind(cor.frame.env.2, cbind(i,
                                                max(abs(mod$residuals)),
                                                median(mod$residuals),
                                                mod$coefficients[1,1],
                                                mod$sigma,
                                                mod$r.squared,
                                                mod$fstatistic[1],
                                                cov(frq, sort.pops.phase2$tmean),
                                                cor(frq, sort.pops.phase2$tmean),
                                                sort(abs(frq-mean(frq)))[21]-sort(abs(frq-mean(frq)))[20]
  )
  )
  setTxtProgressBar(pb, j)
  j<-j+1
}

#fix row/colnames
colnames(cor.frame.env.2)<-c("id","maxresid","medresid","intercept","stdevresid","rsquared","f","cov","cor","outlier")
rownames(cor.frame.env.2) <- c()
#make numeric columns numeric
for (i in 2:10){
  cor.frame.env.2[,i]<-as.numeric(as.character(cor.frame.env.2[,i]))
}

#merge phase1 and phase2 data, retaining only SNPs called in both phases
match.cor.frame<-merge(cor.frame.env, cor.frame.env.2, by='id')
#subtract each column (drop f, cov, and cor due to correlation)
diff.frame<-data.frame(id=match.cor.frame$id,
                       max.resid=abs(match.cor.frame$maxresid.x-match.cor.frame$maxresid.y),
                       med.resid=abs(match.cor.frame$medresid.x-match.cor.frame$medresid.y),
                       int=abs(match.cor.frame$intercept.x-match.cor.frame$intercept.y),
                       stdevresid=abs(match.cor.frame$stdevresid.x-match.cor.frame$stdevresid.y),
                       rsquared=abs(match.cor.frame$rsquared.x-match.cor.frame$rsquared.y),
                       #f=abs(match.cor.frame$f.x-match.cor.frame$f.y),
                       #cov=abs(match.cor.frame$cov.x-match.cor.frame$cov.y),
                       #cor=abs(match.cor.frame$cor.x-match.cor.frame$cor.y),
                       outlier=abs(match.cor.frame$outlier.x-match.cor.frame$outlier.y)
)

#pull out our random subset to hand identify as training dataset
#validate.env.training.set<-diff.frame[sample(nrow(diff.frame),100),]
#validate.env.training.set<-read.csv(file="~/Downloads/bayescenv.temp.validate.env.csv")

par(mfrow=c(3,3))
#plot
sig.pos<-as.vector(validate.env.training.set$id)
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == i,3:16]))
  frq2<-as.numeric(as.vector(phase2.all.freqs[phase2.all.freqs$id == i,3:23]))
  plot(sort.pops.phase1$tmean,frq, main = i, ylim =c(0,1))
  abline(lm(frq~sort.pops.phase1$tmean))
  points(sort.pops.phase2$tmean,frq2, col="red")
  abline(lm(frq2~sort.pops.phase2$tmean), col="red")
  setTxtProgressBar(pb, j)
  j<-j+1
}

#"accept"=1,"reject"=0
vec<-c(1,1,1,1,1,1,0,1,1,
       0,1,1,1,1,1,1,1,1,
       1,1,1,1,1,1,1,1,0,
       1,1,1,1,1,1,1,1,1,
       1,1,1,1,1,1,0,1,1,
       1,1,0,1,1,0,1,1,1,
       1,1,1,1,1,1,0,1,1,
       1,1,1,1,1,1,1,1,0,
       1,1,1,1,1,1,1,1,1,
       1,1,1,1,1,1,1,1,1,
       1,1,1,1,0,1,1,1,0,
       0)

vec[vec == 0]<-"reject"
vec[vec == 1]<-"accept"

#add hand classification to the training set
validate.env.training.set$class<-as.factor(vec)
#write.csv(validate.env.training.set, file="~/Downloads/bayescenv.temp.validate.env.csv")

#train the model (not including max residual which decreased prediction accuracy)
m <- naiveBayes(validate.env.training.set[,2:7], validate.env.training.set[,8])
#check predictions against training data
table(predict(m, validate.env.training.set[,2:7]), validate.env.training.set[,8])
validate.env.training.set$pred.class<-predict(m, validate.env.training.set[,2:7]) #add prediction to df
validate.env.training.set$prob.accept<-predict(m, validate.env.training.set[,2:7], type="raw")[,1] #add prediction probability to df

#do prediction
diff.frame$class<-predict(m, diff.frame[,2:7]) #add prediction to df
diff.frame$prob.accept<-predict(m, diff.frame[,2:7], type="raw")[,1] #add prediction probability to df

#visualize prediction
par(mfrow=c(1,2))
pc.df<-as.data.frame(prcomp(validate.env.training.set[,2:7])$x)
plot(pc.df$PC1,pc.df$PC2, col=validate.env.training.set$class, main = "hand classified")
plot(pc.df$PC1,pc.df$PC2, col=predict(m, validate.env.training.set[,2:7]), main = "ML predicted")

#check prediction
table(diff.frame$class)
par(mfrow=c(1,1))
hist(diff.frame$prob.accept)
#split confidently identified SNPs by class
vetted.bayescenv.temp.env<-diff.frame[diff.frame$prob.accept >.8,]
#add order column from allsigs to visualize this on a manhattan plot
all.sigs.env <- all.sigs[paste(all.sigs$chrom, all.sigs$pos) %in% vetted.bayescenv.temp.env$id,]
all.sigs.env<-data.frame(order=all.sigs.env$order,
                         id=paste(all.sigs.env$chrom,all.sigs.env$pos))
vetted.bayescenv.temp.env<-merge(all.sigs.env,vetted.bayescenv.temp.env,by="id")
write.csv(vetted.bayescenv.temp.env[,1:2], file = "~/Downloads/vetted.bayescenv.temp.env.csv")


#####
####
###
##
#now outliers
#generate matching id's for the outlier dataset in question
sig.pos<-as.vector(cor.frame.outliers$id)
#subset all AF values to only the outlier dataset, to save time rather than searching the entire 3M snp dataset
phase2.all.freqs.bayescenv.temp <- phase2.all.freqs[phase2.all.freqs$id %in% sig.pos,]
#generate list of SNP IDs found in phase2
sig.pos<-phase2.all.freqs.bayescenv.temp$id

#calc phase2 cor.frame
cor.frame.outliers.2<-data.frame(matrix(, nrow=0, ncol=10)) #init df
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase2.all.freqs.bayescenv.temp[phase2.all.freqs.bayescenv.temp$id == i,3:23]))
  mod<-summary(lm(frq ~ sort.pops.phase2$tmean)) #perform linear reg
  # store 9 values describing the correlation to feed into machine learning
  cor.frame.outliers.2<-rbind(cor.frame.outliers.2, cbind(i,
                                                          max(abs(mod$residuals)),
                                                          median(mod$residuals),
                                                          mod$coefficients[1,1],
                                                          mod$sigma,
                                                          mod$r.squared,
                                                          mod$fstatistic[1],
                                                          cov(frq, sort.pops.phase2$tmean),
                                                          cor(frq, sort.pops.phase2$tmean),
                                                          sort(abs(frq-mean(frq)))[21]-sort(abs(frq-mean(frq)))[20]
  )
  )
  setTxtProgressBar(pb, j)
  j<-j+1
}

#fix row/colnames
colnames(cor.frame.outliers.2)<-c("id","maxresid","medresid","intercept","stdevresid","rsquared","f","cov","cor","outlier")
rownames(cor.frame.outliers.2) <- c()
#make numeric columns numeric
for (i in 2:10){
  cor.frame.outliers.2[,i]<-as.numeric(as.character(cor.frame.outliers.2[,i]))
}

#merge phase1 and phase2 data, retaining only SNPs called in both phases
match.cor.frame<-merge(cor.frame.outliers, cor.frame.outliers.2, by='id')
#subtract each column
#subtract each column (drop max.resid, rsquared, and cov due to correlation)
diff.frame<-data.frame(id=match.cor.frame$id,
                       max.resid=abs(match.cor.frame$maxresid.x-match.cor.frame$maxresid.y),
                       med.resid=abs(match.cor.frame$medresid.x-match.cor.frame$medresid.y),
                       int=abs(match.cor.frame$intercept.x-match.cor.frame$intercept.y),
                       stdevresid=abs(match.cor.frame$stdevresid.x-match.cor.frame$stdevresid.y),
                       rsquared=abs(match.cor.frame$rsquared.x-match.cor.frame$rsquared.y),
                       #f=abs(match.cor.frame$f.x-match.cor.frame$f.y),
                       #cov=abs(match.cor.frame$cov.x-match.cor.frame$cov.y),
                       #cor=abs(match.cor.frame$cor.x-match.cor.frame$cor.y),
                       outlier=abs(match.cor.frame$outlier.x-match.cor.frame$outlier.y)
)

#pull out our random subset to hand identify as training dataset
#validate.outliers.training.set<-diff.frame[sample(nrow(diff.frame),100),]
#validate.outliers.training.set<-read.csv(file="~/Downloads/bayescenv.temp.validate.outliers.csv")
#validate.outliers.training.set<-validate.outliers.training.set[,-1]

par(mfrow=c(3,3))
#plot
sig.pos<-as.vector(validate.outliers.training.set$id)
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in sig.pos){
  frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == i,3:16]))
  frq2<-as.numeric(as.vector(phase2.all.freqs[phase2.all.freqs$id == i,3:23]))
  plot(sort.pops.phase1$tmean,frq, main = i, ylim =c(0,1))
  abline(lm(frq~sort.pops.phase1$tmean))
  points(sort.pops.phase2$tmean,frq2, col="red")
  abline(lm(frq2~sort.pops.phase2$tmean), col="red")
  setTxtProgressBar(pb, j)
  j<-j+1
}

#"accept"=1,"reject"=0
vec<-c(0,0,0,0,0,1,0,0,1,
       0,0,0,0,0,0,1,1,0,
       0,1,0,0,0,1,0,0,0,
       1,0,0,1,0,0,0,0,0,
       0,1,1,1,1,0,0,0,0,
       0,0,1,0,0,0,1,0,1,
       0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,1,0,0,
       0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,
       1,1,0,0,1,0,0,0,0,
       0)

vec[vec == 0]<-"reject"
vec[vec == 1]<-"accept"

#add hand classification to the training set
validate.outliers.training.set$class<-as.factor(vec)
#write.csv(validate.outliers.training.set, file="~/Downloads/bayescenv.temp.validate.outliers.csv")

#train the model (not including max residual which decreased prediction accuracy)
m <- naiveBayes(validate.outliers.training.set[,2:7], validate.outliers.training.set[,8])
#check predictions against training data
table(predict(m, validate.outliers.training.set[,2:7]), validate.outliers.training.set[,8])
validate.outliers.training.set$pred.class<-predict(m, validate.outliers.training.set[,2:7]) #add prediction to df
validate.outliers.training.set$prob.accept<-predict(m, validate.outliers.training.set[,2:7], type="raw")[,1] #add prediction probability to df

#do prediction
diff.frame$class<-predict(m, diff.frame[,2:7]) #add prediction to df
diff.frame$prob.accept<-predict(m, diff.frame[,2:7], type="raw")[,1] #add prediction probability to df

#visualize prediction
par(mfrow=c(1,2))
pc.df<-as.data.frame(prcomp(validate.outliers.training.set[,2:7])$x)
plot(pc.df$PC1,pc.df$PC2, col=validate.outliers.training.set$class, main = "hand classified")
plot(pc.df$PC1,pc.df$PC2, col=predict(m, validate.outliers.training.set[,2:7]), main = "ML predicted")

#check prediction
table(diff.frame$class)
par(mfrow=c(1,1))
hist(diff.frame$prob.accept)
#split confidently identified SNPs by class
vetted.bayescenv.temp.outliers<-diff.frame[diff.frame$prob.accept >.8,]
#add order column from allsigs to visualize this on a manhattan plot
all.sigs.out <- all.sigs[paste(all.sigs$chrom, all.sigs$pos) %in% vetted.bayescenv.temp.outliers$id,]
all.sigs.out<-data.frame(order=all.sigs.out$order,
                         id=paste(all.sigs.out$chrom,all.sigs.out$pos))
vetted.bayescenv.temp.outliers<-merge(all.sigs.out,vetted.bayescenv.temp.outliers,by="id")

#pull allele freqs for outliers
phase1.all.freqs.vetted.outliers<- phase1.all.freqs[phase1.all.freqs$id %in% vetted.bayescenv.temp.outliers$id,]

#identify the outlier pop
outlier.pop<-c()
for (i in 1:nrow(phase1.all.freqs.vetted.outliers)){
  outlier.pop[i]<-colnames(phase1.all.freqs.vetted.outliers)[3:16][outlier(as.numeric(phase1.all.freqs.vetted.outliers[i,3:16]), logical = TRUE)]
}
#add to df
vetted.bayescenv.temp.outliers$outlier.pop<-outlier.pop
table(vetted.bayescenv.temp.outliers$outlier.pop)

write.csv(vetted.bayescenv.temp.outliers[,c(1:2,11)], file = "~/Downloads/vetted.bayescenv.temp.outliers.csv")

#









