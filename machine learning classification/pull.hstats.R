#pull out hstats files



temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_AOM.txt.gz", temp)
gzfile(temp, 'rt')
AO.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
AO.2R$order<-AO.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_AOM.txt.gz", temp)
gzfile(temp, 'rt')
AO.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
AO.2L$order<-AO.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_AOM.txt.gz", temp)
gzfile(temp, 'rt')
AO.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
AO.3R$order<-AO.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_AOM.txt.gz", temp)
gzfile(temp, 'rt')
AO.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
AO.3L$order<-AO.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_AOM.txt.gz", temp)
gzfile(temp, 'rt')
AO.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
AO.X$order<-AO.X$start+61545105+49364325+53200684+41963435
#
AO<-as.data.frame(rbind(AO.2R,AO.2L,AO.3R,AO.3L,AO.X))
write.csv(AO, file="~/Downloads/AO.hstats.csv")



temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_GAS.txt.gz", temp)
gzfile(temp, 'rt')
GA.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GA.2R$order<-GA.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_GAS.txt.gz", temp)
gzfile(temp, 'rt')
GA.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GA.2L$order<-GA.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_GAS.txt.gz", temp)
gzfile(temp, 'rt')
GA.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GA.3R$order<-GA.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_GAS.txt.gz", temp)
gzfile(temp, 'rt')
GA.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GA.3L$order<-GA.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_GAS.txt.gz", temp)
gzfile(temp, 'rt')
GA.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GA.X$order<-GA.X$start+61545105+49364325+53200684+41963435
#
GA<-as.data.frame(rbind(GA.2R,GA.2L,GA.3R,GA.3L,GA.X))
write.csv(GA, file="~/Downloads/GA.hstats.csv")


temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_GWA.txt.gz", temp)
gzfile(temp, 'rt')
GW.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GW.2R$order<-GW.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_GWA.txt.gz", temp)
gzfile(temp, 'rt')
GW.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GW.2L$order<-GW.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_GWA.txt.gz", temp)
gzfile(temp, 'rt')
GW.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GW.3R$order<-GW.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_GWA.txt.gz", temp)
gzfile(temp, 'rt')
GW.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GW.3L$order<-GW.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_GWA.txt.gz", temp)
gzfile(temp, 'rt')
GW.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
GW.X$order<-GW.X$start+61545105+49364325+53200684+41963435
#
GW<-as.data.frame(rbind(GW.2R,GW.2L,GW.3R,GW.3L,GW.X))
write.csv(GW, file="~/Downloads/GW.hstats.csv")


temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_KES.txt.gz", temp)
gzfile(temp, 'rt')
KE.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
KE.2R$order<-KE.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_KES.txt.gz", temp)
gzfile(temp, 'rt')
KE.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
KE.2L$order<-KE.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_KES.txt.gz", temp)
gzfile(temp, 'rt')
KE.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
KE.3R$order<-KE.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_KES.txt.gz", temp)
gzfile(temp, 'rt')
KE.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
KE.3L$order<-KE.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_KES.txt.gz", temp)
gzfile(temp, 'rt')
KE.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
KE.X$order<-KE.X$start+61545105+49364325+53200684+41963435
#
KE<-as.data.frame(rbind(KE.2R,KE.2L,KE.3R,KE.3L,KE.X))
write.csv(KE, file="~/Downloads/KE.hstats.csv")


temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_UGS.txt.gz", temp)
gzfile(temp, 'rt')
UG.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
UG.2R$order<-UG.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_UGS.txt.gz", temp)
gzfile(temp, 'rt')
UG.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
UG.2L$order<-UG.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_UGS.txt.gz", temp)
gzfile(temp, 'rt')
UG.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
UG.3R$order<-UG.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_UGS.txt.gz", temp)
gzfile(temp, 'rt')
UG.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
UG.3L$order<-UG.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_UGS.txt.gz", temp)
gzfile(temp, 'rt')
UG.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
UG.X$order<-UG.X$start+61545105+49364325+53200684+41963435
#
UG<-as.data.frame(rbind(UG.2R,UG.2L,UG.3R,UG.3L,UG.X))
write.csv(UG, file="~/Downloads/UG.hstats.csv")


temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_BFM.txt.gz", temp)
gzfile(temp, 'rt')
BFM.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFM.2R$order<-BFM.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_BFM.txt.gz", temp)
gzfile(temp, 'rt')
BFM.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFM.2L$order<-BFM.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_BFM.txt.gz", temp)
gzfile(temp, 'rt')
BFM.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFM.3R$order<-BFM.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_BFM.txt.gz", temp)
gzfile(temp, 'rt')
BFM.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFM.3L$order<-BFM.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_BFM.txt.gz", temp)
gzfile(temp, 'rt')
BFM.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFM.X$order<-BFM.X$start+61545105+49364325+53200684+41963435
#
BFM<-as.data.frame(rbind(BFM.2R,BFM.2L,BFM.3R,BFM.3L,BFM.X))
write.csv(BFM, file="~/Downloads/BFM.hstats.csv")


temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2R_BFS.txt.gz", temp)
gzfile(temp, 'rt')
BFS.2R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFS.2R$order<-BFS.2R$start
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_2L_BFS.txt.gz", temp)
gzfile(temp, 'rt')
BFS.2L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFS.2L$order<-BFS.2L$start+61545105
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3R_BFS.txt.gz", temp)
gzfile(temp, 'rt')
BFS.3R <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFS.3R$order<-BFS.3R$start+61545105+49364325
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_3L_BFS.txt.gz", temp)
gzfile(temp, 'rt')
BFS.3L <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFS.3L$order<-BFS.3L$start+61545105+49364325+53200684
#
temp <- tempfile()
download.file("ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/selection/hstats/hstats_table_X_BFS.txt.gz", temp)
gzfile(temp, 'rt')
BFS.X <- read.table(temp, sep = "\t",header = T)[,c(2,3,6)]
unlink(temp)
BFS.X$order<-BFS.X$start+61545105+49364325+53200684+41963435
#
BFS<-as.data.frame(rbind(BFS.2R,BFS.2L,BFS.3R,BFS.3L,BFS.X))
write.csv(BFS, file="~/Downloads/BFS.hstats.csv")







