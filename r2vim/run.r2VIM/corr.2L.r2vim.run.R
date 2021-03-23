library(randomForest)
library(Pomona)
library(ranger)
library(parallel)

##read in adjusted df
reg_data_corrected<-read.csv("2L.corr.corrected.csv")
p<-ncol(reg_data_corrected)-2

#read in position info
pos<-read.csv("2L.chrom.no.missing.pos.csv")

#execute 10 random forests with 1000 trees each, and the optimized mtry value: via this wrapper function for #ranger
r2vim<-var.sel.r2vim(x = reg_data_corrected[,3:ncol(reg_data_corrected)], y=reg_data_corrected$response,
                     no.runs = 10, factor = 1, ntree = 1000, mtry.prop = .1,
                     method = "ranger", type = "regression")

#calculate correlation between runs to ensure ntree=1000 was high enough to achieve convergence
runs<-as.data.frame(cbind(r2vim$info$vim.run.1,r2vim$info$vim.run.2,r2vim$info$vim.run.3,r2vim$info$vim.run.4,r2vim$info$vim.run.5,
                          r2vim$info$vim.run.6,r2vim$info$vim.run.7,r2vim$info$vim.run.8,r2vim$info$vim.run.9,r2vim$info$vim.run.10))
#calc correlation
cor(runs)
jpeg(file="corrs.2L.corr.jpeg")
hist(cor(runs), col="blue")
dev.off()
plot(runs[,1],runs[,2]) #compare run 1 and 2 visually

#make results df for manhattan plot
df<-data.frame(chrom=pos$chrom,
               snp=pos$bp,
               min.rel.vim=r2vim$info$rel.vim.min)

#plot
plot(x=1:nrow(df),y=r2vim$info$rel.vim.min)
abline(h=1, col="red")

#save rel.vim.min values
write.csv(df, "2L.corr.rel.vim.min.csv", row.names=FALSE, quote=FALSE)
