
#combine all significance values into a single df
#already converted lfmm p vals to q vals and combined runs into single df
#already combined bayescenv q vals into single df
#already combined rel.vim.mins into single df
#now stick them all together

lfmm<-read.csv("~/Downloads/lfmm.qvalues.csv")
bayescenv<-read.csv("~/Desktop/anoph.phase2/bayescenv.results/baye.qvalues.output.csv")
r2vim<-read.csv("~/Desktop/anoph.phase2/r2vim/all.vars.rel.vim.mins.csv")
rep.corrs<-read.csv("~/Desktop/anoph.phase2/repeated.corrs.csv")
table(paste(lfmm$Chromosome,lfmm$Position) == paste(bayescenv$chrom, bayescenv$pos))
table(paste(lfmm$Chromosome,lfmm$Position) == paste(r2vim$chrom, r2vim$snp))
table(paste(lfmm$Chromosome,lfmm$Position) == paste(rep.corrs$chrom, rep.corrs$pos))
rep.corrs<-rep.corrs[match(paste(lfmm$Chromosome,lfmm$Position),paste(rep.corrs$chrom, rep.corrs$pos)),]
table(paste(lfmm$Chromosome,lfmm$Position) == paste(rep.corrs$chrom, rep.corrs$pos))
head(lfmm)
head(bayescenv)
head(r2vim)
table(bayescenv == 0)
bayescenv[bayescenv == 0]<-.00002

all.sigs<-data.frame(Chromosome=lfmm$Chromosome,
                     Position=lfmm$Position,
                     lfmm.prec=lfmm$p_value_p_annual,
                     lfmm.hum=lfmm$p_value_h_annual,
                     lfmm.temp=lfmm$p_value_t_mean,
                     lfmm.corr=lfmm$p_value_autocorrelated,
                     lfmm.rand=lfmm$p_value_random,
                     baye.prec=bayescenv$prec.q,
                     baye.hum=bayescenv$hum.q,
                     baye.temp=bayescenv$temp.q,
                     baye.corr=bayescenv$corr.q,
                     baye.rand=bayescenv$rand.q,
                     vim.prec=r2vim$vim.prec,
                     vim.hum=r2vim$vim.hum,
                     vim.temp=r2vim$vim.temp,
                     vim.corr=r2vim$vim.corr,
                     vim.rand=r2vim$vim.rand,
                     cors.prec=rep.corrs$prec,
                     cors.hum=rep.corrs$hum,
                     cors.temp=rep.corrs$temp,
                     cors.corr=rep.corrs$corr,
                     cors.rand=rep.corrs$rand)


write.csv(all.sigs, file="~/Downloads/all.sigs.csv", row.names = F, quote = F)
