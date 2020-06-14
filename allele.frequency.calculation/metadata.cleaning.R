library(dplyr)
library(ggplot2)

#bring in dataframe with locality and environmental data for each individual
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")

#bring in phase2 metadata
phase2.alldata<-read.csv(file = "~/Downloads/phase2.metadata.csv")

#awesome dplyr function that keeps all the occurrences in the first dataset not found in the second dataset
phase2.uniques<-anti_join(phase2.alldata, phase1.alldata, by = 'ox_code')


#use dplyr to make a table with all unique lat longs and their frequency phase2
site.sampling.phase2<-phase2.uniques %>% group_by(latitude, longitude) %>% summarize(count=n())

#use dplyr to make a table with all unique lat longs and their frequency phase2
site.sampling.phase1<-phase1.alldata %>% group_by(latitude, longitude) %>% summarize(count=n())
#use dplyr to make a table with all unique lat longs and their frequency phase2
site.sampling.phase2.alldata<-phase2.alldata %>% group_by(latitude, longitude) %>% summarize(count=n())

#define africa map
afr<-map_data("world", xlim=c(-20,60),ylim=c(-35,40))  

#map phase1
ggplot()+
  scale_y_continuous(limits = c(-35,35)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = site.sampling.phase1, aes(x = longitude, y = latitude), color = "red", size = 3.5, alpha =.5) +
  ggtitle("new samples phase1")
#map phase1 by sample size
ggplot()+
  scale_y_continuous(limits = c(-35,35)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = site.sampling.phase1, aes(x = longitude, y = latitude), color = "red", size = .15*site.sampling.phase1$count, alpha =.5) +
  ggtitle("new samples phase1")


#map phase2
ggplot()+
  scale_y_continuous(limits = c(-35,35)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = phase2.uniques, aes(x = longitude, y = latitude), color = "red", size = 3.5, alpha =.1) +
  ggtitle("new samples phase2")
#map phase2 size by num samples
ggplot()+
  scale_y_continuous(limits = c(-35,35)) +
  scale_x_continuous(limits = c(-25,55)) +
  geom_polygon(data = afr, aes(x=long, y = lat, group = group))+
  geom_point(data = site.sampling.phase2, aes(x = longitude, y = latitude), color = "red", size = .15*site.sampling.phase2$count, alpha =.5) +
  ggtitle("new samples phase2")


#questions:
#Is it okay to modify points from mayotte? What should the concensus point be? Most frequent lat/long?
#Is it okay to combine (-3.86200, 39.74500) with (-3.63500, 39.85800) ?
#Is it okay to drop (6.09449, -0.26093)

#shift all mayotte points to (-12.852147, 45.10389)
phase2.uniques$latitude[phase2.uniques$latitude == -12.990662 | phase2.uniques$latitude == -12.857 |
                        phase2.uniques$latitude == -12.796525 | phase2.uniques$latitude == -12.778704 |
                        phase2.uniques$latitude == -12.737813 | phase2.uniques$latitude == -12.70271] <- -12.852147
phase2.uniques$longitude[phase2.uniques$longitude > 45] <- 45.10389

#shift (-3.63500, 39.85800) to (-3.86200, 39.74500)
phase2.uniques$latitude[phase2.uniques$latitude == -3.63500] <- -3.86200
phase2.uniques$longitude[phase2.uniques$longitude == 39.85800] <- 39.74500

#drop (6.09449, -0.26093)
phase2.uniques$latitude == 6.09449
phase2.uniques<-phase2.uniques[c(1:45,47:377),]


#write file of individuals to keep to filter vcf file by
write.table(phase2.uniques$ox_code, file="~/Downloads/phase2.uniques.txt", row.names = F, quote = F, col.names = F)

#write csv of locality info for included phase2 dataset
write.csv(phase2.uniques[,c(1,15,16)], file="~/Downloads/phase2.unique.sampling.locs.csv", row.names = F)

#next:
#create plots of hum vs precip vs temp for each method to check correlation btwn variables
#create outlier manhattan plots
#create allele frequency plots with phase 2 overlaid







