
#make sampling figure for the paper

#read in phase1 data
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
#read in unique phase2 points
phase2.alldata<-read.csv(file = "~/Downloads/phase2.unique.sampling.locs.csv")

#subset matching columns
ph1<-phase1.alldata[,colnames(phase2.alldata)]

#use dplyr to make a table with all unique lat longs and their frequency phase2
ph1.site.sampling<-ph1 %>% group_by(latitude, longitude) %>% summarize(count=n())

#use dplyr to make a table with all unique lat longs and their frequency phase2
ph2.site.sampling<-phase2.alldata %>% group_by(latitude, longitude) %>% summarize(count=n())

#add phase column
ph2.site.sampling$phase<-rep(2, times=nrow(ph2.site.sampling))
#add phase column
ph1.site.sampling$phase<-rep(1, times=nrow(ph1.site.sampling))

#combine in tidy format
all.phase<-as.data.frame(rbind(ph1.site.sampling,ph2.site.sampling))
all.phase$phase<-as.factor(all.phase$phase)


##define africa map
afr<-map_data("world", xlim=c(-20,60),ylim=c(-35,40))  

#map phase1
ggplot()+
  scale_y_continuous(limits = c(-35,37)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = all.phase, aes(x = longitude, y = latitude, col = phase), size = 3.5, alpha =.5) +
  ggtitle("new samples phase1")+
  theme_bw()


#map phase1 by sample size
ggplot()+
  scale_y_continuous(limits = c(-35,37)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = all.phase, aes(x = longitude, y = latitude, col=phase), size = .15*all.phase$count, alpha =.6) +
  theme_bw()

#map phase1 by sample size
ggplot()+
  scale_y_continuous(limits = c(-35,37)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = all.phase, aes(x = longitude, y = latitude, col=phase, size=count), alpha =.6) +
  theme_bw()

#map phase1 by sample size
ggplot()+
  scale_y_continuous(limits = c(-35,37)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = all.phase, aes(x = longitude, y = latitude, col=phase, size=count), alpha =.75) +
  theme_classic()

ggsave(filename = "~/Downloads/anoph.sampling.map.pdf",width = 168, height = 119.5, units="mm")
ggsave(filename = "~/Downloads/anoph.sampling.map.jpeg",width = 168, height = 119.5, units="mm")








