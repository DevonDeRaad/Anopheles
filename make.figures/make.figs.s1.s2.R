library(ggplot2)
#make Fig S1
x<-c(253,99,102,272,171,144,830,1211,178)
y<-c(1,0,0,7,7,3,9,9,0)
df<-data.frame(genes=x,go.terms=y)
ggplot(df, aes(x=genes, y=go.terms)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  theme_classic()+
  labs(x="unique genes used for GO term enrichment testing",
       y="significantly enriched GO terms")+
  scale_y_continuous(breaks=seq(0, 10, by = 2))+
  scale_x_continuous(breaks=seq(0, 1200, by = 300))
ggsave("~/Desktop/anoph.phase2/fig.s1.pdf", width = 4, height = 4, units = "in")
ggsave("~/Desktop/anoph.phase2/fig.s1.png", width = 4, height = 4, units = "in")

#make Fig S2
#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c("latitude","longitude","pannual","tmean","hannual")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase1 <-metapops[order(metapops$latitude),]

ggplot(sort.pops.phase1, aes(x=pannual, y=hannual)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  theme_classic()+
  labs(x="Precipitation",
       y="Humidity")

ggsave("~/Desktop/anoph.phase2/fig.s2.pdf", width = 3, height = 3, units = "in")
ggsave("~/Desktop/anoph.phase2/fig.s2.png", width = 3, height = 3, units = "in")
