#day6 NMDS
#load multiple pkgs
lapply(c("vegan", "cluster"), require, character.only = TRUE)

#source
source('biostats.R')

speabu <- read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <- read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1) 
# Let’s transform the data before diving into the analysis. 
# We will do this because the species abundance dataset is highly 
# skewed and contains some rather large values which are valid but highly influential.
# Remember, the log of zero is undefined so we’ll add 1 to each value in our data set.

speabu.log<-log(speabu+1)

#Let’s perform metaMDS using our log-transformed dataset.
spe.nmds<-metaMDS(speabu.log, distance='bray', k=2, autotransform=FALSE, trymax=100)
spe.nmds


names(spe.nmds)

spe.nmds2<-metaMDS(speabu.log, distance='bray', k=3, autotransform=FALSE, trymax=100)
# Examining a scree plot of stress versus the number of dimensions can help you make this decision. 
# To perform this type:
nmds.scree(speabu.log, distance='bray', k=10, autotransform=FALSE, trymax=20)

#a Monte Carlo randomization test of the final stress value can be conducted as follows. Note that this will take a couple minutes to complete.
nmds.monte(speabu.log, distance='bray', k=3, autotransform=FALSE, trymax=20)
#This will return the permuted stress values (and histogram) and calculated p-value.

stressplot(spe.nmds2)

plot(spe.nmds,type='n') 
text(spe.nmds,labels=row.names(speabu))

#We can make the symbol size proportional to log abundance. Try typing,
plot(spe.nmds,type='n') 
points(spe.nmds,cex=speabu.log$ROSYDACE)

#Scores
spe.nmds$points

#Calculate the loadings (i.e., variable weights)

#A permutation test is used to assess statistical significance, rather than using the F distribution.
vec.sp<-envfit(spe.nmds$points, speabu.log, perm=1000)

vec.sp

#
#plot these loadings on the ordination plot.
par()
ordiplot(spe.nmds, choices = c(1, 2), type="text", display="sites", xlab="Axis 1", ylab="Axis 2")
plot(vec.sp, p.max=.01, col="blue")

#plotting , by product


Env_variables<- colnames(envdata)

#default multipanel plot
# 3 rows and 4 columns - mf - multipanel plot
par(mfcol=c(3,4),  oma=c(5,3,1,1),mar=c(1,3,0,0))
for (i in 1:10) {
  plot(x=spe.nmds$points[,'MDS1'], y=envdata[,c(i)], type="p",  
       xaxs="i",yaxs="i")        #remove the first two columns which contain data
  abline(lm(envdata[,c(i)] ~ spe.nmds$points[,'MDS1']),col='red')
  abline(1,1)
  box(col="gray80") 
  
  
  
  #side -3  top, line = -1
  mtext(side=3,line=-1,text=Env_variables[i],cex=0.5) #add area label
  
}