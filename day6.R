#day6
# FISH 560: Applied Multivariate Statistics for Ecologists
# T opics
# Principal coordinate analysis (PCoA) R Packages: vegan
# R Source: biostats

#load multiple pkgs
lapply(c("vegan"), require, character.only = TRUE)
#source
source('biostats.R')


speabu <- read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
# Let’s transform the data before diving into the analysis.
# We will do this because the species abundance dataset is 
# highly skewed and contains some rather large values which are valid 
# but highly influential. Remember, the log of zero is undefined so
# we’ll add 1 to each value in our data set.
speabu.log<-log(speabu+1)

# Calculating the distance matrix
# The first step in PCoA is is the selection of a distance matrix to 
# describe similarities among objects according to the descriptors


speabu.d<-vegdist(speabu.log, "bray")
#Let’s perform PCoA:
spe.pcoa<-cmdscale(speabu.d, eig=TRUE, add=T)

#The principal scores are contained in spe.pcoa$points 
#and the eigenvalues are contained in 
spe.pcoa$eig

#Let’s calculate this for the first five PCs.
spe.pcoa$eig/sum(spe.pcoa$eig)*100

#eigenvalues to expectations according to the broken stick model.
plot(spe.pcoa$eig[1:35]/sum(spe.pcoa$eig)*100,type="b",lwd=2,col="blue",
     xlab= "Principal Component from PCoA", 
     ylab="% variation explained", main="% variation explained by PCoA (blue) vs. random expectation (red)")
lines(bstick(35)*100,type="b",lwd=2,col="red")


# view the ordination plot:
ordiplot(spe.pcoa, choices = c(1, 2), type="text",
         display="sites", xlab="PC 1 (15%)", ylab="PC 2 (11%)")

#Calculate the PC loadings (i.e., variable weights)
vec.sp<-envfit(spe.pcoa$points, k=45, speabu.log, perm=1000)
#Another way of doing this is to extract the first two PCs from our PCoA 
#object using the scores() function. To do this, type:
vec.sp<-envfit(scores(spe.pcoa), speabu.log, perm=1000)
vec.sp

names(vec.sp) 
names(vec.sp$vectors)

par()
ordiplot(spe.pcoa, choices = c(1, 2), type="text", display="sites", xlab="PC 1 (15%)", ylab="PC 2 (11%)")
plot(vec.sp, p.max=.01, col="blue")
